#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
import argparse
import re
import sys

DEFAULT_NEW_FILE = Path("details_new.txt")
DEFAULT_OLD_FILE = Path("details.txt")
DEFAULT_OUT_FILE = Path("details_new.txt")


try:
    from ftlib import Compound2 as _FTLIB_Compound2
except Exception as exc:
    _FTLIB_IMPORT_ERROR = exc

    class _FTLIB_Compound2:
        pass

Compound2 = _FTLIB_Compound2
sys.modules["__main__"].Compound2 = Compound2


def expand_ads_token(tok):
    tok = tok.strip()
    if not tok:
        return []

    m = re.match(r"^(\d+)(.+)$", tok)
    if m:
        return [m.group(2)] * int(m.group(1))
    return [tok]


def canonical_ref(ref_parts):
    if not ref_parts:
        return Counter()

    # Spezialfall: "A 4" am Zeilenende
    if len(ref_parts) == 2 and ref_parts[1].isdigit():
        return Counter({ref_parts[0]: int(ref_parts[1])})

    text = " ".join(ref_parts).strip()
    if not text:
        return Counter()

    tokens = []
    for part in text.split("+"):
        tokens.extend(expand_ads_token(part))

    return Counter(tokens)


def ref_to_string(ref_counter):
    if not ref_counter:
        return ""

    parts = []
    for ads in sorted(ref_counter):
        n = ref_counter[ads]
        parts.append(ads if n == 1 else f"{n}{ads}")

    return "+".join(parts)


def normalize_ads_for_match(x):
    return str(x).strip().rstrip("/").split("/")[-1]


def ads_suffix_match(short_ads, full_ads):
    """
    Erste Spalte: nur exakter Name oder exaktes Ende mit '+'-Grenze.
    Keine chemische Kanonisierung.
    """
    short_ads = normalize_ads_for_match(short_ads)
    full_ads = normalize_ads_for_match(full_ads)

    if short_ads == full_ads:
        return True

    return full_ads.endswith("+" + short_ads)


def ads_prefix_before_suffix(short_ads, full_ads):
    """
    Return the prefix that was omitted when short_ads is used as suffix of full_ads.

    Examples:
      short C, full H2O+H2O+C -> H2O+H2O
      short CH, full CH4+4CH4+C+CH -> CH4+4CH4+C
      short 4CH4_ts, full CH4+4CH4_ts -> CH4
    """
    short_ads = normalize_ads_for_match(short_ads)
    full_ads = normalize_ads_for_match(full_ads)

    if short_ads == full_ads:
        return ""

    suffix = "+" + short_ads
    if full_ads.endswith(suffix):
        return full_ads[: -len(suffix)]

    return None


def strip_whole_ts_suffix_for_sequence(text):
    """
    For sequence/ref-context checks only, treat final _ts as a marker for the
    whole adsorbate name.

    This fixes:
      CH4+4CH4_ts -> CH4+4CH4
      4CH4_ts     -> 4CH4

    The original names are still kept for exact first-column output matching.
    """
    text = str(text).strip()
    if text.endswith("_ts"):
        return text[:-3]
    return text


def ts_suffix_alias_allowed(old_item, new_item):
    """
    Special TS alias:
      old first col: 4CH4_ts
      new first col: CH4+4CH4_ts

    The normal suffix guard can reject this because the omitted prefix CH4 is
    not necessarily part of the manual reference. For TS rows, however the
    short name is commonly the actual TS tail while the generated name contains
    the full context.
    """
    old_ads = normalize_ads_for_match(old_item["ads"])
    new_ads = normalize_ads_for_match(new_item["ads"])

    if not old_ads.endswith("_ts"):
        return False

    omitted = ads_prefix_before_suffix(old_ads, new_ads)
    if omitted in (None, ""):
        return False

    # Do not allow hiding another TS in the omitted context.
    if "_ts" in omitted:
        return False

    return True


def counter_subset(a, b):
    return all(a[k] <= b[k] for k in a)


def ads_sequence(text):
    """
    Ordered expanded sequence for prefix checks.

    This is intentionally stricter than Counter subset:
      H2O+2H2O+C -> [H2O, H2O, H2O, C]
      H2O+3H2O+C -> [H2O, H2O, H2O, H2O, C]

    Therefore H2O+2H2O+C is NOT a prefix of H2O+3H2O+C.

    For TS names, a final _ts is treated as marker of the whole string:
      CH4+4CH4_ts -> [CH4, CH4, CH4, CH4, CH4]
      4CH4_ts     -> [CH4, CH4, CH4, CH4]
    """
    text = strip_whole_ts_suffix_for_sequence(text)
    if not text:
        return []

    parts = text.split()

    # Accept reference shorthand: "A 4"
    if len(parts) == 2 and parts[1].isdigit():
        return [parts[0]] * int(parts[1])

    seq = []
    for part in text.split("+"):
        seq.extend(expand_ads_token(part))
    return seq


def omitted_prefix_is_covered_by_reference(omitted_prefix, old_item, new_item):
    """
    Guard for first-column short names.

    A short first-column name is valid only if:

      1) the omitted prefix is the ordered beginning of the MANUAL reference, and
      2) any remaining part of the manual reference is explicitly present at the
         beginning of the manual first-column adsorbate name.

    This is the missing distinction:

    Valid:
      manual first col: 3H2O+C+O
      manual ref:       H2O+3H2O
      auto first col:   H2O+3H2O+C+O

      omitted prefix = H2O
      remaining manual ref after omitted prefix = 3H2O
      manual first col starts with 3H2O -> valid.

    Invalid:
      manual first col: CH4
      manual ref:       H2O+4H2O
      auto first col:   H2O+3H2O+CH4

      omitted prefix = H2O+3H2O
      remaining manual ref after omitted prefix = H2O
      manual first col does not start with H2O -> invalid.

    Invalid:
      manual first col: C+C2H4
      manual ref:       H2O+3H2O+C
      auto first col:   H2O+2H2O+C+C+C2H4

      omitted prefix is not even the ordered beginning of manual ref -> invalid.
    """
    if omitted_prefix in (None, ""):
        return omitted_prefix == ""

    omitted_seq = ads_sequence(omitted_prefix)
    manual_ref_seq = ads_sequence(old_item.get("ref_raw", ""))
    manual_ads_seq = ads_sequence(old_item.get("ads", ""))

    if not omitted_seq:
        return True

    if len(omitted_seq) > len(manual_ref_seq):
        return False

    if manual_ref_seq[: len(omitted_seq)] != omitted_seq:
        return False

    remaining_manual_ref = manual_ref_seq[len(omitted_seq):]

    if not remaining_manual_ref:
        return True

    if len(remaining_manual_ref) > len(manual_ads_seq):
        return False

    return manual_ads_seq[: len(remaining_manual_ref)] == remaining_manual_ref


def suffix_alias_allowed(old_item, new_item):
    """
    True if old first-column name may be used as a short name for new first-column name.

    This keeps exact matches untouched, but prevents cases like:
      old: CH                  ref CH4+4CH4
      new: CH4+4CH4+C+CH       ref CH4+4CH4

    Here the short name CH would hide the extra C, so it is not allowed.
    """
    old_ads = old_item["ads"]
    new_ads = new_item["ads"]

    if old_ads == new_ads:
        return True

    if not ads_suffix_match(old_ads, new_ads):
        return False

    # Special case for TS aliases such as:
    #   4CH4_ts -> CH4+4CH4_ts
    # Keep the strict guarded logic for non-TS names.
    if ts_suffix_alias_allowed(old_item, new_item):
        return True

    omitted = ads_prefix_before_suffix(old_ads, new_ads)
    return omitted_prefix_is_covered_by_reference(omitted, old_item, new_item)


def first_column_name_compatible(a, b):
    """
    Symmetrische Namens-Kompatibilität für Ketten-Zwischenschritte.

    This is deliberately looser than suffix_alias_allowed because the chain
    search later decides by atom-set union and final reference. Do not use this
    for the main old/new alias assignment.
    """
    return ads_suffix_match(a, b) or ads_suffix_match(b, a)


def parse_line(line):
    parts = line.split()
    if len(parts) < 2:
        return None

    ads = parts[0]
    nums_raw = parts[1]
    ref_parts = parts[2:]
    ref_raw = " ".join(ref_parts).strip()

    try:
        nums = tuple(sorted(int(x) for x in nums_raw.split(",") if x.strip()))
    except ValueError:
        return None

    return ads, nums, canonical_ref(ref_parts), ref_raw


def read_details(path):
    entries = {}
    lines = []

    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")

    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue

        parsed = parse_line(stripped)
        if parsed is None:
            print(f"WARNING: cannot parse {path}:{lineno}: {stripped}")
            lines.append(stripped)
            continue

        ads, nums, ref_counter, ref_raw = parsed

        if ads in entries:
            print(f"WARNING: duplicate adsorbate in {path}: {ads}")

        entries[ads] = {
            "ads": ads,
            "nums": nums,
            "nums_set": set(nums),
            "ref": ref_counter,
            "ref_raw": ref_raw,
            "line": stripped,
            "lineno": lineno,
            "source": str(path),
        }

        lines.append(stripped)

    return entries, lines


def entry_key(item):
    return (item["source"], item["lineno"], item["ads"], item["line"])


def build_adsorbate_alias_mapping(old_entries, new_entries):
    """
    Baut eine vorsichtige 1:1-Zuordnung:

      new_ads -> old_ads

    Für die erste Spalte gilt absichtlich KEINE chemische Kanonisierung.
    Erlaubt ist:
      1) exakter Name
      2) eindeutiges exaktes String-Ende mit '+'-Grenze,
         ABER nur wenn der weggelassene Prefix durch die Referenz abgedeckt ist.

    Dadurch wird verhindert:
      CH  -> CH4+4CH4+C+CH  bei gleicher Referenz CH4+4CH4

    Denn dann müsste manual entweder C+CH heißen oder CH mit der längeren
    Referenz CH4+4CH4+C.
    """
    old_ads = list(old_entries.keys())
    new_ads = list(new_entries.keys())

    new_to_old = {}
    old_to_new = {}

    def assign(new, old):
        new_to_old[new] = old
        old_to_new[old] = new

    # 1) Exakte Matches zuerst. Alte Logik unverändert.
    for new in new_ads:
        if new in old_entries and new not in new_to_old and new not in old_to_new:
            assign(new, new)

    # 2) Eindeutige, aber referenz-geschützte Suffix-Matches iterativ.
    progress = True
    while progress:
        progress = False

        remaining_old = [old for old in old_ads if old not in old_to_new]
        remaining_new = [new for new in new_ads if new not in new_to_old]

        # Spezifische Namen zuerst, z.B. 3H2O+C+O vor C+O.
        remaining_old.sort(key=len, reverse=True)

        for old in remaining_old:
            if old in old_to_new:
                continue

            old_item = old_entries[old]
            cands = []
            for new in remaining_new:
                new_item = new_entries[new]
                if suffix_alias_allowed(old_item, new_item):
                    cands.append(new)

            if len(cands) == 1:
                assign(cands[0], old)
                progress = True

        if progress:
            continue

        remaining_old = [old for old in old_ads if old not in old_to_new]
        remaining_new = [new for new in new_ads if new not in new_to_old]

        remaining_new.sort(key=len, reverse=True)

        for new in remaining_new:
            if new in new_to_old:
                continue

            new_item = new_entries[new]
            cands = []
            for old in remaining_old:
                old_item = old_entries[old]
                if suffix_alias_allowed(old_item, new_item):
                    cands.append(old)

            if len(cands) == 1:
                assign(new, cands[0])
                progress = True

    ambiguous = []
    blocked = []

    remaining_old = [old for old in old_ads if old not in old_to_new]
    remaining_new = [new for new in new_ads if new not in new_to_old]

    for old in remaining_old:
        raw_cands = [
            new for new in remaining_new
            if ads_suffix_match(old, new)
        ]
        allowed_cands = [
            new for new in raw_cands
            if suffix_alias_allowed(old_entries[old], new_entries[new])
        ]

        if len(allowed_cands) > 1:
            ambiguous.append((old, allowed_cands))

        for new in raw_cands:
            if new not in allowed_cands:
                blocked.append((old, new))

    return new_to_old, ambiguous, blocked


def candidate_entries_for_name(query_ads, new_entries, old_entries):
    """
    Liefert ALLE Kandidaten, nicht nur eindeutige.

    Das ist der zentrale Fix:
    Bei query H2O+3H2O+C dürfen sowohl 3H2O+C als auch C Kandidaten sein.
    Entschieden wird danach über Referenz + Atom-Summen-Test.
    """
    query_ads = normalize_ads_for_match(query_ads)
    candidates = []

    # generierte Details bevorzugt, danach manuelle
    for source_name, entries in [("new", new_entries), ("old", old_entries)]:
        # exact first
        if query_ads in entries:
            candidates.append((source_name + ":exact", entries[query_ads]))

        for ads, item in entries.items():
            if ads == query_ads:
                continue
            if first_column_name_compatible(query_ads, ads):
                candidates.append((source_name + ":suffix", item))

    # deduplicate exact/suffix duplicates
    out = []
    seen = set()
    for how, item in candidates:
        key = entry_key(item)
        if key not in seen:
            out.append((how, item))
            seen.add(key)

    return out


def try_reference_chain_match(old_item, new_item, new_entries, old_entries, max_depth=6):
    """
    Akzeptiert Fälle, in denen generated eine längere Referenz benutzt als manual.

    Es gibt zwei erlaubte Fälle:

    A) Die Zwischenstruktur erklärt zusätzliche Atome der manuellen Zeile.
       Dann wird union(new atoms + intermediate atoms) == manual atoms geprüft.

    B) Die eigentlichen Adsorbatatome passen bereits exakt.
       Dann dürfen Zwischenstrukturen nur als Referenzkontext-Brücke dienen.
       Ihre Atomnummern werden NICHT zur Ziel-Atomliste addiert.

       Genau dieser Fall:
         manual: OH                69,74  CH4+3CH4
         auto:   CH4+3CH4+C+OH     69,74  CH4+3CH4+C
         ref context CH4+3CH4+C kann in details z.B. als C 72 CH4+3CH4 stehen.
         Atom 72 darf hier ignoriert werden, weil OH selbst schon stimmt.
    """
    target_ref = old_item["ref"]
    target_nums = set(old_item["nums"])

    start_nums = set(new_item["nums"])
    if not start_nums.issubset(target_nums):
        return False, [f'{new_item["ads"]}: start nums not subset of manual nums']

    def suffix_after_prefix(prefix, full):
        """
        Return the first-column tail of full after prefix.

        Example:
          prefix CH4+3CH4+C
          full   CH4+3CH4+C+OH
          -> OH
        """
        prefix = normalize_ads_for_match(prefix)
        full = normalize_ads_for_match(full)

        if not prefix or not full:
            return None

        if full == prefix:
            return ""

        needle = prefix + "+"
        if full.startswith(needle):
            return full[len(needle):]

        return None

    def new_tail_is_manual_tail():
        """
        The auto line may use a longer reference and therefore write only the
        final adsorbate tail.

        Example:
          auto   CH4+3CH4+C+OH 69,74 CH4+3CH4+C
          manual 3CH4+C+OH     64,69,74 CH4+3CH4

        The actual adsorbate tail is OH. If auto nums are a subset of manual
        nums and the reference context is valid, the extra manual atom can be a
        differently numbered reference-extension atom.
        """
        tail = suffix_after_prefix(new_item["ref_raw"], new_item["ads"])

        if tail is None or tail == "":
            return False

        old_ads = normalize_ads_for_match(old_item["ads"])

        return old_ads == tail or old_ads.endswith("+" + tail)

    def ref_context_reaches_target_by_name(query_ref_raw, chain, seen, depth):
        """
        Follow only the reference context, ignoring all candidate atom numbers.

        This function answers:
          Does query_ref_raw correspond to some details entry whose own reference
          eventually reaches target_ref?

        It is deliberately only used when the actual adsorbate atom list already
        matches the manual atom list.
        """
        if depth >= max_depth:
            return False, chain + ["context max-depth"]

        candidates = candidate_entries_for_name(query_ref_raw, new_entries, old_entries)

        if not candidates:
            return False, chain + [f"{query_ref_raw} (context no candidates)"]

        failed = []

        # Prefer candidates whose own reference is already the manual reference.
        def context_sort_key(pair):
            _how, item = pair
            return (0 if item["ref"] == target_ref else 1, len(item["ads"]))

        for how, cand in sorted(candidates, key=context_sort_key):
            key = entry_key(cand)
            if key in seen:
                continue

            step = (
                f"{cand['ads']}({how}; reference-context; "
                f"nums ignored={','.join(map(str, cand['nums']))}; ref={cand['ref_raw']})"
            )

            if cand["ref"] == target_ref:
                return True, chain + [step]

            next_ref = cand["ref_raw"]
            if not next_ref:
                failed.append(step + " -> empty-ref")
                continue

            ok, result = ref_context_reaches_target_by_name(
                next_ref,
                chain + [step],
                seen | {key},
                depth + 1,
            )
            if ok:
                return True, result

            failed.extend(result[-3:])

        return False, chain + failed[:8]

    # Critical early bridge:
    # If the auto atom list is already a subset of the manual atom list and
    # both lines describe the same final tail, we only need to prove that the
    # auto reference context reaches the manual reference. Do NOT add the
    # reference-context atoms, because their indices may differ between manual
    # and auto references.
    #
    # This fixes:
    #   auto   CH4+3CH4+C+OH 69,74    CH4+3CH4+C
    #   manual 3CH4+C+OH     64,69,74 CH4+3CH4
    # where C may be index 72 in the auto reference context but index 64 in the
    # manual shorter-reference line. OH itself still matches as 69,74.
    if start_nums.issubset(target_nums) and new_tail_is_manual_tail() and new_item["ref"] != target_ref:
        ok, context_chain = ref_context_reaches_target_by_name(
            new_item["ref_raw"],
            [
                f"{new_item['ads']}(new-start; nums={','.join(map(str, new_item['nums']))}; ref={new_item['ref_raw']}; tail-only)"
            ],
            {entry_key(new_item)},
            0,
        )
        if ok:
            return True, context_chain
        # If context bridge fails, continue with the old union logic below.
        early_failed = context_chain
    else:
        early_failed = []

    def dfs(current, union_nums, chain, seen, depth):
        if current["ref"] == target_ref:
            if union_nums == target_nums:
                return True, chain
            return False, chain + [
                "ref reached but nums differ "
                f"union={','.join(map(str, sorted(union_nums)))} "
                f"target={','.join(map(str, sorted(target_nums)))}"
            ]

        if depth >= max_depth:
            return False, chain + ["max-depth"]

        next_query = current["ref_raw"]
        if not next_query:
            return False, chain + ["empty-ref"]

        candidates = candidate_entries_for_name(next_query, new_entries, old_entries)

        if not candidates:
            return False, chain + [f"{next_query} (no candidates)"]

        failed_notes = list(early_failed[-4:]) if depth == 0 else []

        def cand_sort_key(pair):
            _how, item = pair
            extra = set(item["nums"]) - target_nums
            adds = len(set(item["nums"]) - union_nums)
            return (len(extra), -adds, len(item["ads"]))

        candidates = sorted(candidates, key=cand_sort_key)

        for how, cand in candidates:
            key = entry_key(cand)
            if key in seen:
                continue

            cand_nums = set(cand["nums"])
            new_union = union_nums | cand_nums

            if not new_union.issubset(target_nums):
                # Secondary bridge fallback:
                # If the real adsorbate atom list is already complete, ignore
                # reference-context candidate atoms and only validate its reference chain.
                if union_nums == target_nums:
                    if cand["ref"] == target_ref:
                        return True, chain + [
                            f"{cand['ads']}({how}; reference-context; nums ignored={','.join(map(str, cand['nums']))}; ref={cand['ref_raw']})"
                        ]

                    if cand["ref_raw"]:
                        ok, context_chain = ref_context_reaches_target_by_name(
                            cand["ref_raw"],
                            chain + [
                                f"{cand['ads']}({how}; reference-context; nums ignored={','.join(map(str, cand['nums']))}; ref={cand['ref_raw']})"
                            ],
                            seen | {key},
                            depth + 1,
                        )
                        if ok:
                            return True, context_chain

                failed_notes.append(
                    f"{cand['ads']}({how}) adds outside target: "
                    f"{','.join(map(str, sorted(new_union - target_nums)))}"
                )
                continue

            ok, result_chain = dfs(
                cand,
                new_union,
                chain + [
                    f"{cand['ads']}({how}; nums={','.join(map(str, cand['nums']))}; ref={cand['ref_raw']})"
                ],
                seen | {key},
                depth + 1,
            )

            if ok:
                return True, result_chain

            failed_notes.extend(result_chain[-3:])

        return False, chain + failed_notes[:10]

    return dfs(
        new_item,
        start_nums,
        [
            f"{new_item['ads']}(new-start; nums={','.join(map(str, new_item['nums']))}; ref={new_item['ref_raw']})"
        ],
        {entry_key(new_item)},
        0,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-file", type=Path, default=DEFAULT_NEW_FILE)
    parser.add_argument("--old-file", type=Path, default=DEFAULT_OLD_FILE)
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUT_FILE)
    parser.add_argument("--debug-chain", action="store_true")
    args = parser.parse_args()

    old_entries, old_lines = read_details(args.old_file)
    new_entries, _ = read_details(args.new_file)

    new_to_old, ambiguous, blocked_aliases = build_adsorbate_alias_mapping(
        old_entries,
        new_entries,
    )

    for old, cands in ambiguous:
        print(
            "WARNING: ambiguous short adsorbate name "
            f"{old!r} could match: {', '.join(cands)}"
        )

    for old, new in blocked_aliases:
        if args.debug_chain:
            print(
                "DEBUG: blocked short-name alias because omitted prefix / remaining manual reference does not match first-column shortcut: "
                f"{old!r} -> {new!r}; "
                f"manual ref={ref_to_string(old_entries[old]['ref'])}; "
                f"auto ref={ref_to_string(new_entries[new]['ref'])}"
            )

    out_lines = list(old_lines)

    n_ok = 0
    n_ok_chain = 0
    n_warn = 0
    n_add = 0
    n_alias = 0

    for ads, item in new_entries.items():
        old_ads = new_to_old.get(ads)

        if old_ads is not None:
            old_item = old_entries[old_ads]

            alias_note = ""
            if old_ads != ads:
                alias_note = f" (manual name {old_ads!r})"
                n_alias += 1

            if old_item["nums"] == item["nums"] and old_item["ref"] == item["ref"]:
                print(f"OK: {ads}{alias_note} numbers/ref match")
                n_ok += 1
            else:
                chain_ok, chain = try_reference_chain_match(
                    old_item=old_item,
                    new_item=item,
                    new_entries=new_entries,
                    old_entries=old_entries,
                )

                if chain_ok:
                    print(f"OK_CHAIN: {ads}{alias_note} matches manual shorter reference")
                    if args.debug_chain:
                        print(f"  chain: {' -> '.join(chain)}")
                    n_ok_chain += 1
                else:
                    print(f"WARNING: mismatch for {ads}{alias_note}")
                    print(f"  original nums: {','.join(map(str, old_item['nums']))}")
                    print(f"  new      nums: {','.join(map(str, item['nums']))}")
                    print(f"  original ref : {ref_to_string(old_item['ref'])}")
                    print(f"  new      ref : {ref_to_string(item['ref'])}")
                    print(f"  tried chain  : {' -> '.join(chain)}")
                    print(f"  original line: {old_item['line']}")
                    print(f"  new      line: {item['line']}")
                    n_warn += 1
        else:
            print(f"ADD: {ads}")
            out_lines.append(item["line"])
            n_add += 1

    args.output.write_text("\n".join(out_lines) + "\n")

    print()
    print(f"Wrote {args.output}")
    print(f"OK:       {n_ok}")
    print(f"OK chain: {n_ok_chain}")
    print(f"Aliases:  {n_alias}")
    print(f"Warnings: {n_warn}")
    print(f"Added:    {n_add}")


if __name__ == "__main__":
    main()
