#!/usr/bin/env python3
"""
list_single_adsorbate_and_ts_dnames_v2.py

Input:
  1) outputs pickle: list of output objects
  2) compounds pickle: dict of compound objects, key = compound.nickname

Ziel:
  Kompakte Python-Liste von dnames, nach Referenz gruppiert:

      dnames = (
          ['Ni_111_C+' + b for b in ['H', 'O', 'CH-H_ts']] +
          ['Ni_111_CO+' + b for b in ['H', 'O']]
      )

Filter:
  - Alle TS-Compounds werden aufgenommen, auch wenn adsorbate "+" enthält.
  - Nicht-TS-Compounds werden nur aufgenommen, wenn der Teil zusätzlich zur Referenz
    genau ein einzelnes Adsorbat ist.
  - comp_ref None / "None" / "" zählt als leere Referenz:
      Prefix = metal_facet_
      Suffix = adsorbate

Wichtig:
  - output.dname selbst hat kein "_ts".
  - Für diese Ausgabeliste wird bei TS am Ende "_ts" angehängt.
"""

from __future__ import annotations

import argparse
import os
import pickle
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

try:
    from ftlib import Compound2
except Exception as exc:
    linear_sum_assignment = None
    _FTLIB_IMPORT_ERROR = exc

def load_pickle(path: Path) -> Any:
    try:
        with path.open("rb") as f:
            return pickle.load(f)
    except Exception as pickle_error:
        try:
            import dill  # type: ignore
        except Exception:
            raise RuntimeError(
                f"Could not unpickle {path} with pickle, and dill is not available.\n"
                f"pickle error:\n{pickle_error}"
            ) from pickle_error

        with path.open("rb") as f:
            return dill.load(f)


def req_attr(obj: Any, attr: str, context: str) -> Any:
    if not hasattr(obj, attr):
        raise AttributeError(f"{context} has no required attribute {attr!r}: {obj!r}")
    return getattr(obj, attr)


def truthy(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, int):
        return x != 0
    if isinstance(x, float):
        return x != 0.0
    if isinstance(x, str):
        return x.strip().lower() in {"true", "t", "yes", "y", "1"}
    return bool(x)


def is_ts(comp: Any) -> bool:
    return truthy(getattr(comp, "ts", False))


def basename_no_slash(x: Any) -> str:
    return os.path.basename(str(x).rstrip("/"))


def is_empty_ref(x: Any) -> bool:
    if x is None:
        return True
    s = str(x).strip()
    return s == "" or s.lower() == "none"


def normalize_name(x: Any) -> str:
    if is_empty_ref(x):
        return ""
    s = basename_no_slash(x).strip()
    if s.endswith("_ts"):
        s = s[:-3]
    return s


def split_adsorbate_tokens(x: Any) -> list[str]:
    s = normalize_name(x)
    if not s:
        return []
    return [tok.strip() for tok in s.split("+") if tok.strip()]


def token_counter(x: Any) -> Counter[str]:
    return Counter(split_adsorbate_tokens(x))


def additional_tokens(adsorbate: Any, comp_ref: Any) -> list[str]:
    """
    Entfernt die Referenz-Tokens einmalig aus adsorbate.
    Beispiel:
      adsorbate='C+H', comp_ref='C' -> ['H']
      adsorbate='4CH4+CH3+H', comp_ref='4CH4' -> ['CH3', 'H']
    """
    ads_tokens = split_adsorbate_tokens(adsorbate)
    ref_counts = token_counter(comp_ref)

    rest: list[str] = []
    for tok in ads_tokens:
        if ref_counts[tok] > 0:
            ref_counts[tok] -= 1
        else:
            rest.append(tok)

    return rest


def compound_nickname(key: Any, comp: Any) -> str:
    return str(getattr(comp, "nickname", key))


def metal_facet_prefix(key: Any, comp: Any) -> str:
    metal = str(req_attr(comp, "metal", f"compound {key!r}"))
    facet = str(req_attr(comp, "facet", f"compound {key!r}"))
    return f"{metal}_{facet}_"


def expected_output_name(key: Any, comp: Any) -> str:
    # output.dname has no "_ts"
    return metal_facet_prefix(key, comp) + str(req_attr(comp, "adsorbate", f"compound {key!r}"))


class OutputResolver:
    def __init__(self, outputs: list[Any]):
        self.by_dname: dict[str, str] = {}

        for out in outputs:
            if not hasattr(out, "dname"):
                continue

            dname = str(getattr(out, "dname"))
            self.by_dname[dname] = dname
            self.by_dname[basename_no_slash(dname)] = dname

    def canonical(self, dname: str) -> str | None:
        if dname in self.by_dname:
            return self.by_dname[dname]
        return self.by_dname.get(basename_no_slash(dname))


def output_exists_for_compound(key: Any, comp: Any, outputs: OutputResolver) -> str | None:
    """
    Prüft, ob der Compound in outputs vorhanden ist.
    Rückgabe ist actual output.dname, ohne "_ts".
    """
    candidates = [
        str(key),
        compound_nickname(key, comp),
        basename_no_slash(key),
        basename_no_slash(compound_nickname(key, comp)),
        expected_output_name(key, comp),
        basename_no_slash(expected_output_name(key, comp)),
    ]

    seen: set[str] = set()
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)

        out_dname = outputs.canonical(cand)
        if out_dname is not None:
            return out_dname

    return None


def keep_compound(comp: Any) -> bool:
    """
    TS immer behalten.
    Nicht-TS nur, wenn der Zusatz zur Referenz genau ein Token ist.
    """
    if is_ts(comp):
        return True

    ads = req_attr(comp, "adsorbate", "compound")
    ref = getattr(comp, "comp_ref", None)
    rest = additional_tokens(ads, ref)

    return len(rest) == 1 and "+" not in rest[0]


def reference_prefix_and_suffix(key: Any, comp: Any) -> tuple[str, str]:
    """
    Gruppierung nach Referenz.

    Beispiele:
      metal=Ni, facet=111, comp_ref=None, ads=C
        -> prefix 'Ni_111_', suffix 'C'

      comp_ref=C, ads=C+H
        -> prefix 'Ni_111_C+', suffix 'H'

      comp_ref=4CH4, ads=4CH4+CH3+H, ts=True
        -> prefix 'Ni_111_4CH4+', suffix 'CH3+H_ts'
    """
    base_prefix = metal_facet_prefix(key, comp)
    ads = normalize_name(req_attr(comp, "adsorbate", f"compound {key!r}"))
    ref = getattr(comp, "comp_ref", None)
    ref_norm = normalize_name(ref)

    rest = additional_tokens(ads, ref)

    if ref_norm:
        prefix = base_prefix + ref_norm + "+"
        suffix = "+".join(rest) if rest else ads
    else:
        prefix = base_prefix
        suffix = ads

    if is_ts(comp) and not suffix.endswith("_ts"):
        suffix += "_ts"

    return prefix, suffix


def format_python_expression(groups: dict[str, list[str]], varname: str) -> str:
    parts: list[str] = []

    for prefix, suffixes in sorted(groups.items(), key=lambda kv: kv[0]):
        suffixes = sorted(set(suffixes))
        parts.append(f"[{prefix!r} + b for b in {suffixes!r}]")

    if not parts:
        expr = "[]"
    else:
        expr = " +\n    ".join(parts)

    return f"{varname} = (\n    {expr}\n)"


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--outputs_pickle",
        type=Path,
        default=Path("structures.pkl"),
        help="Pickle file containing the output objects. Default: structures.pkl",
    )

    parser.add_argument(
        "--compounds_pickle",
        type=Path,
        default=Path("compounds.pkl"),
        help="Pickle file containing the compound dict. Default: compounds.pkl",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional: write expression to file instead of only stdout.",
    )

    parser.add_argument(
        "--varname",
        default="dnames",
        help="Variable name in printed Python expression. Default: dnames",
    )

    parser.add_argument(
        "--debug-skips",
        action="store_true",
        help="Print compounds without matching output.dname to stderr.",
    )

    args = parser.parse_args()

    outputs = load_pickle(args.outputs_pickle)
    compounds = load_pickle(args.compounds_pickle)

    if not isinstance(outputs, list):
        raise TypeError(f"{args.outputs_pickle} must contain a list, got {type(outputs).__name__}")
    if not isinstance(compounds, dict):
        raise TypeError(f"{args.compounds_pickle} must contain a dict, got {type(compounds).__name__}")

    output_resolver = OutputResolver(outputs)

    groups: dict[str, list[str]] = defaultdict(list)

    skipped_no_output = 0
    skipped_filter = 0

    for key, comp in compounds.items():
        if not keep_compound(comp):
            skipped_filter += 1
            continue

        out_dname = output_exists_for_compound(key, comp, output_resolver)
        if out_dname is None:
            skipped_no_output += 1
            if args.debug_skips:
                print(f"SKIP no output.dname match for compound key={key!r}", file=__import__("sys").stderr)
            continue

        prefix, suffix = reference_prefix_and_suffix(key, comp)
        groups[prefix].append(suffix)

    expression = format_python_expression(groups, args.varname)

    print(expression)
    print()
    print(f"# selected: {sum(len(set(v)) for v in groups.values())}")
    print(f"# reference groups: {len(groups)}")
    print(f"# skipped by filter: {skipped_filter}")
    print(f"# skipped without output.dname match: {skipped_no_output}")

    if args.output is not None:
        args.output.write_text(expression + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
