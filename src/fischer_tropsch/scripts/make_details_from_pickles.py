#!/usr/bin/env python3
"""
print_adsorbate_rows_with_refs_write_exact_v11_safe_output.py

Zurück zur funktionierenden Tabellen-Idee, nur minimal erweitert.

Was dieses Skript macht:

1) CompoundData baut beim __init__ eine echte Liste:
       self.rows

2) Jede Row enthält:
       row.key              # nur für Summenformel/Matching
       row.key_formula
       row.adsorbate        # erste Spalte in details_new.txt
       row.ref
       row.dname
       row.atoms
       row.adsorbate_dict
       row.ref_row
       row.ref_atoms
       row.ref_adsorbate_dict

3) --print-table druckt wieder die einfache Tabelle:
       key key_formula adsorbate ref dname adsorbate_dict

4) --test-exact druckt nur einen einfachen Test:
       Gibt es in row.adsorbate_dict mindestens ein Adsorbat
       mit der Summenformel row.key_formula?

Schreibt jetzt details_new.txt für exakte Fälle und details_new.log nur mit --write-log als nachvollziehbares Log.
Noch keine Zerlegung in mehrere Adsorbate.
"""

import argparse
import os
import pickle
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


# Pickle compatibility for old __main__.Compound2 pickles.
try:
    from ftlib import Compound2 as _FTLIB_Compound2
except Exception:
    class _FTLIB_Compound2:
        pass

Compound2 = _FTLIB_Compound2
sys.modules["__main__"].Compound2 = Compound2


class Row:
    def __init__(self, key, key_formula, adsorbate, ref, dname, atoms, adsorbate_dict,
                 metal, facet, source):
        # key is only used for formula matching / exact selection.
        self.key = key
        self.key_formula = key_formula

        # adsorbate is the name written as first column to details_new.txt.
        # For normal compounds this is compound.adsorbate.
        # For missing reference rows this is derived from the output dname.
        self.adsorbate = adsorbate

        self.ref = ref
        self.dname = dname
        self.atoms = atoms
        self.adsorbate_dict = adsorbate_dict
        self.metal = metal
        self.facet = facet
        self.source = source

        # Wird nach dem Tabellenaufbau gesetzt.
        self.ref_row = None
        self.ref_atoms = None
        self.ref_adsorbate_dict = {}


class CompoundData:
    COV_RAD = {
        "H": 0.31,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "S": 1.05,
        "P": 1.07,
        "F": 0.57,
        "Cl": 1.02,
        "Br": 1.20,
        "I": 1.39,
    }

    def __init__(self, outputs, compounds, surface_elements=None, bond_scale=1.35, ts_bond_scale=1.75):
        self.outputs = list(outputs.values()) if isinstance(outputs, dict) else list(outputs)
        self.compounds = compounds
        self.surface_override = surface_elements
        self.bond_scale = bond_scale
        self.ts_bond_scale = ts_bond_scale

        self.outputs_by_dname = self.build_outputs_by_dname()
        self.compounds_by_dname = self.build_compounds_by_dname()
        # Für "existiert in compounds?" zählen ALLE Compounds.
        # Auch comp_ref=None Compounds verhindern also, dass eine missing-ref-output
        # Row nachgebaut wird. In self.rows kommen sie trotzdem nicht.
        self.compound_keys = set(self.bname(k) for k in self.compounds)

        # Wichtig: das ist die ganze Tabelle.
        self.rows = []
        self.rows_by_dname = {}

        self.build_rows()
        self.link_reference_rows()

    @classmethod
    def from_files(cls, outputs_pickle, compounds_pickle, surface_elements=None, bond_scale=1.35, ts_bond_scale=1.75):
        outputs = cls.load_pickle(outputs_pickle)
        compounds = cls.load_pickle(compounds_pickle)
        if not isinstance(compounds, dict):
            raise TypeError("compounds pickle must contain a dict, got %s" % type(compounds).__name__)
        return cls(
            outputs,
            compounds,
            surface_elements=surface_elements,
            bond_scale=bond_scale,
            ts_bond_scale=ts_bond_scale,
        )

    @staticmethod
    def load_pickle(path):
        try:
            with open(path, "rb") as f:
                return pickle.load(f)
        except Exception:
            import dill
            with open(path, "rb") as f:
                return dill.load(f)

    @staticmethod
    def opt(obj, attr, default=None):
        return getattr(obj, attr, default)

    @staticmethod
    def bname(x):
        if x is None:
            return ""
        s = os.path.basename(str(x).strip().rstrip("/"))
        return "" if s.lower() == "none" else s

    @staticmethod
    def is_atoms(obj):
        return obj is not None and hasattr(obj, "get_chemical_symbols") and hasattr(obj, "get_positions")

    @staticmethod
    def truthy(x):
        if isinstance(x, bool):
            return x
        if isinstance(x, str):
            return x.strip().lower() in ("1", "true", "yes", "y", "t")
        return bool(x)

    @classmethod
    def strip_ts(cls, name):
        s = cls.bname(name)
        return s[:-3] if s.endswith("_ts") else s

    @staticmethod
    def unique(items):
        out = []
        seen = set()
        for x in items:
            if x and x not in seen:
                seen.add(x)
                out.append(x)
        return out

    @classmethod
    def dname_variants(cls, dname):
        return cls.unique([dname, dname.replace("-", "+"), dname.replace("+", "-")])

    # ----- names -----

    def metal(self, comp):
        return str(self.opt(comp, "metal", "")).strip()

    def facet(self, comp):
        return str(self.opt(comp, "facet", "")).strip()

    def adsorbate(self, comp):
        return self.bname(self.opt(comp, "adsorbate", ""))

    def compound_ref_name(self, comp):
        return self.bname(self.opt(comp, "comp_ref", None))

    def compound_has_reference(self, comp):
        return bool(self.compound_ref_name(comp))

    def is_co_pool_name(self, name):
        name = self.strip_ts(name).replace("-", "+").strip()
        if not name:
            return False

        total = 0
        for token in name.split("+"):
            token = token.strip()
            if not token:
                continue

            m = re.match(r"^(\d*)CO$", token)
            if not m:
                return False

            n = int(m.group(1)) if m.group(1) else 1
            total += n

        return total >= 2

    def expected_dname(self, comp):
        return "%s_%s_%s" % (self.metal(comp), self.facet(comp), self.strip_ts(self.adsorbate(comp)))

    def expected_ref_dname(self, metal, facet, ref):
        return "%s_%s_%s" % (metal, facet, self.strip_ts(ref))

    def expected_ref_dname_from_comp(self, parent_comp, ref):
        return self.expected_ref_dname(self.metal(parent_comp), self.facet(parent_comp), ref)

    def output_dname_candidates(self, comp):
        ads = self.adsorbate(comp)
        candidates = self.dname_variants(self.expected_dname(comp))

        if ads.endswith("_ts"):
            base = ads[:-3]
            candidates += self.dname_variants("%s_%s_%s" % (self.metal(comp), self.facet(comp), base))

        return self.unique(candidates)

    def output_dname_candidates_for_ref(self, parent_comp, ref):
        ref_name = self.strip_ts(ref)
        candidates = self.dname_variants(
            "%s_%s_%s" % (self.metal(parent_comp), self.facet(parent_comp), ref_name)
        )

        raw = self.bname(ref)
        if raw.endswith("_ts"):
            candidates += self.dname_variants(
                "%s_%s_%s" % (self.metal(parent_comp), self.facet(parent_comp), raw[:-3])
            )

        return self.unique(candidates)

    def ads_part_from_dname(self, dname, comp):
        prefix = "%s_%s_" % (self.metal(comp), self.facet(comp))
        if dname.startswith(prefix):
            return dname[len(prefix):]
        return dname

    def key_ref_from_output_dname(self, dname, comp):
        ads_part = self.ads_part_from_dname(dname, comp)
        ads_part = self.strip_ts(ads_part).replace("-", "+")

        parts = [p for p in ads_part.split("+") if p]
        if not parts:
            return self.bname(dname), ""

        key = parts[-1]
        ref = "+".join(parts[:-1])
        return key, ref

    def adsorbate_from_output_dname(self, dname, comp):
        ads_part = self.ads_part_from_dname(dname, comp)
        return self.strip_ts(ads_part).replace("-", "+")

    # ----- indices -----

    def build_outputs_by_dname(self):
        out = {}
        for obj in self.outputs:
            dname = self.bname(self.opt(obj, "dname", ""))
            if dname:
                out[dname] = obj
        return out

    def build_compounds_by_dname(self):
        out = defaultdict(list)
        for key, comp in self.compounds.items():
            # Für reine Existenzchecks zählen ALLE Compounds.
            # Ob sie in self.rows kommen, entscheidet build_rows().
            for dname in self.dname_variants(self.expected_dname(comp)):
                out[dname].append((key, comp))
        return dict(out)

    def find_output(self, candidates):
        for dname in candidates:
            if dname in self.outputs_by_dname:
                return self.outputs_by_dname[dname], dname
        return None, candidates[0] if candidates else ""

    def find_output_for_compound(self, comp):
        return self.find_output(self.output_dname_candidates(comp))

    def reference_exists_in_compounds(self, ref, parent_comp):
        ref_name = self.bname(ref)
        if not ref_name:
            return True

        if ref_name in self.compound_keys:
            return True

        for dname in self.dname_variants(self.expected_ref_dname_from_comp(parent_comp, ref_name)):
            if dname in self.compounds_by_dname:
                return True

        return False

    # ----- formula -----

    def formula_from_counter(self, c):
        order = []
        for el in ("C", "H"):
            if c.get(el, 0):
                order.append(el)
        for el in sorted(c):
            if el not in order and c[el]:
                order.append(el)

        return "".join(el if c[el] == 1 else "%s%d" % (el, c[el]) for el in order)

    def parse_formula_token(self, token):
        token = self.strip_ts(token)
        if "_" in token:
            token = token.split("_", 1)[0]

        if not token:
            return None

        parts = re.findall(r"([A-Z][a-z]?)(\d*)", token)
        if not parts or "".join(el + n for el, n in parts) != token:
            return None

        c = Counter()
        for el, n in parts:
            c[el] += int(n) if n else 1
        return c

    def key_formula(self, key):
        name = self.strip_ts(key).replace("-", "+")
        total = Counter()

        for token in name.split("+"):
            token = token.strip()
            if not token:
                continue

            m = re.match(r"^(\d+)([A-Z].*)$", token)
            if m:
                mult = int(m.group(1))
                token = m.group(2)
            else:
                mult = 1

            c = self.parse_formula_token(token)
            if c is None:
                continue

            for el, n in c.items():
                total[el] += mult * n

        return self.formula_from_counter(total) if total else ""

    def expanded_key_unit_formulas(self, key):
        """
        Zerlegt row.key in einzelne Adsorbat-Einheiten.

        Beispiele:
          A+B   -> [A, B]
          3A    -> [A, A, A]
          A+2B  -> [A, B, B]

        Zurückgegeben werden Summenformeln der Einheiten.
        """
        name = self.strip_ts(key).replace("-", "+")
        formulas = []

        for token in name.split("+"):
            token = token.strip()
            if not token:
                continue

            m = re.match(r"^(\d+)([A-Z].*)$", token)
            if m:
                mult = int(m.group(1))
                token = m.group(2)
            else:
                mult = 1

            c = self.parse_formula_token(token)
            if c is None:
                continue

            formula = self.formula_from_counter(c)
            for _ in range(mult):
                formulas.append(formula)

        return formulas

    # ----- atoms / adsorbate dict -----

    def choose_atoms(self, comp, output_obj):
        if output_obj is not None and self.is_atoms(self.opt(output_obj, "transition_state", None)):
            return output_obj.transition_state
        if output_obj is not None and self.is_atoms(self.opt(output_obj, "atoms", None)):
            return output_obj.atoms
        if self.is_atoms(self.opt(comp, "atoms", None)):
            return comp.atoms
        return None

    def choose_output_atoms(self, output_obj):
        if output_obj is not None and self.is_atoms(self.opt(output_obj, "transition_state", None)):
            return output_obj.transition_state
        if output_obj is not None and self.is_atoms(self.opt(output_obj, "atoms", None)):
            return output_obj.atoms
        return None

    def add_ts_suffix_if_needed(self, name, is_ts):
        if is_ts and name and not name.endswith("_ts"):
            return name + "_ts"
        return name

    def row_adsorbate_from_compound(self, comp, output_obj):
        name = self.adsorbate(comp)
        is_ts = (
            self.truthy(self.opt(comp, "ts", False))
            or (output_obj is not None and self.is_atoms(self.opt(output_obj, "transition_state", None)))
            or name.endswith("_ts")
        )
        return self.add_ts_suffix_if_needed(name, is_ts)

    def row_adsorbate_from_output(self, dname, comp, output_obj):
        name = self.adsorbate_from_output_dname(dname, comp)
        is_ts = output_obj is not None and self.is_atoms(self.opt(output_obj, "transition_state", None))
        return self.add_ts_suffix_if_needed(name, is_ts)

    def surface_elements_for_comp(self, comp):
        if self.surface_override:
            return set(x.strip() for x in self.surface_override.split(",") if x.strip())
        return set(re.findall(r"[A-Z][a-z]?", self.metal(comp)))

    def surface_elements_for_metal(self, metal):
        if self.surface_override:
            return set(x.strip() for x in self.surface_override.split(",") if x.strip())
        return set(re.findall(r"[A-Z][a-z]?", metal))

    def covrad(self, symbol):
        return self.COV_RAD.get(symbol, 0.80)

    def mic_distance(self, atoms, i, j):
        pos = atoms.get_positions()
        vec = pos[j] - pos[i]
        try:
            from ase.geometry import find_mic
            vec, _ = find_mic(vec, cell=atoms.get_cell().array, pbc=atoms.get_pbc())
        except Exception:
            pass
        return float((vec * vec).sum() ** 0.5)

    def connected_components(self, atoms, surface_elements, is_ts=False):
        symbols = atoms.get_chemical_symbols()
        idxs = [i for i, s in enumerate(symbols) if s not in surface_elements]

        adj = {i: set() for i in idxs}

        scale = self.ts_bond_scale if is_ts else self.bond_scale

        for n, i in enumerate(idxs):
            for j in idxs[n + 1:]:
                cutoff = scale * (self.covrad(symbols[i]) + self.covrad(symbols[j]))
                if self.mic_distance(atoms, i, j) <= cutoff:
                    adj[i].add(j)
                    adj[j].add(i)

        comps = []
        seen = set()

        for i in idxs:
            if i in seen:
                continue

            stack = [i]
            seen.add(i)
            comp = []

            while stack:
                k = stack.pop()
                comp.append(k)
                for nb in adj[k]:
                    if nb not in seen:
                        seen.add(nb)
                        stack.append(nb)

            comps.append(sorted(comp))

        return comps

    def molecule_formula(self, atoms, indices):
        symbols = atoms.get_chemical_symbols()
        return self.formula_from_counter(Counter(symbols[i] for i in indices))

    def adsorbate_dict_from_atoms(self, atoms, metal, is_ts=False):
        if atoms is None:
            return {}

        surface = self.surface_elements_for_metal(metal)
        out = defaultdict(list)

        for indices in self.connected_components(atoms, surface, is_ts=is_ts):
            out[self.molecule_formula(atoms, indices)].append(indices)

        return dict(sorted(out.items()))

    # ----- rows -----

    def add_row(self, row):
        self.rows.append(row)

        # Index nur zum Finden der Referenz.
        # Die Tabelle selbst bleibt self.rows.
        if row.dname and row.dname not in self.rows_by_dname:
            self.rows_by_dname[row.dname] = row

    def make_row(self, key, adsorbate, ref, dname, atoms, metal, facet, source):
        is_ts = self.bname(adsorbate).endswith("_ts")
        return Row(
            key=key,
            key_formula=self.key_formula(key),
            adsorbate=adsorbate,
            ref=ref,
            dname=dname,
            atoms=atoms,
            adsorbate_dict=self.adsorbate_dict_from_atoms(atoms, metal, is_ts=is_ts),
            metal=metal,
            facet=facet,
            source=source,
        )

    def compound_row(self, key, comp):
        output_obj, dname = self.find_output_for_compound(comp)
        atoms = self.choose_atoms(comp, output_obj)
        ref = self.bname(self.opt(comp, "comp_ref", None))

        return self.make_row(
            key=self.bname(key),
            adsorbate=self.row_adsorbate_from_compound(comp, output_obj),
            ref=ref,
            dname=dname,
            atoms=atoms,
            metal=self.metal(comp),
            facet=self.facet(comp),
            source="compound",
        )

    def missing_reference_row(self, ref, parent_comp):
        output_obj, dname = self.find_output(self.output_dname_candidates_for_ref(parent_comp, ref))
        key, ref_part = self.key_ref_from_output_dname(dname, parent_comp)
        adsorbate = self.row_adsorbate_from_output(dname, parent_comp, output_obj)
        atoms = self.choose_output_atoms(output_obj)

        return self.make_row(
            key=key,
            adsorbate=adsorbate,
            ref=ref_part,
            dname=dname,
            atoms=atoms,
            metal=self.metal(parent_comp),
            facet=self.facet(parent_comp),
            source="missing-ref-output",
        )

    def build_rows(self):
        # 1) Normale Compound-Zeilen.
        # Compounds mit comp_ref None/leer werden komplett ausgeschlossen.
        for key, comp in self.compounds.items():
            if not self.compound_has_reference(comp):
                continue

            if self.is_co_pool_name(key):
                continue

            row = self.compound_row(key, comp)
            self.add_row(row)

        # 2) Fehlende Referenzen als gleichberechtigte Zeilen.
        printed_missing_refs = set()

        for _key, comp in self.compounds.items():
            if not self.compound_has_reference(comp):
                continue

            ref = self.opt(comp, "comp_ref", None)
            ref_name = self.bname(ref)
            if not ref_name:
                continue

            if self.reference_exists_in_compounds(ref, comp):
                continue

            marker = (self.metal(comp), self.facet(comp), ref_name)
            if marker in printed_missing_refs:
                continue
            printed_missing_refs.add(marker)

            row = self.missing_reference_row(ref, comp)
            self.add_row(row)

    def link_reference_rows(self):
        for row in self.rows:
            row.ref_row = None
            row.ref_atoms = None
            row.ref_adsorbate_dict = {}

            if not row.ref:
                continue

            candidates = self.dname_variants(self.expected_ref_dname(row.metal, row.facet, row.ref))
            for dname in candidates:
                ref_row = self.rows_by_dname.get(dname)
                if ref_row is not None:
                    row.ref_row = ref_row
                    row.ref_atoms = ref_row.atoms
                    row.ref_adsorbate_dict = ref_row.adsorbate_dict
                    break

    # ----- print/debug -----

    def print_table(self):
        for row in self.rows:
            print(row.key, row.key_formula, row.adsorbate, row.ref, row.dname, row.adsorbate_dict)

    def print_exact_formula_test(self):
        for i, row in enumerate(self.rows):
            matches = row.adsorbate_dict.get(row.key_formula, []) if row.key_formula else []
            ref_ok = row.ref_row is not None

            print(
                "ROW",
                i,
                "key=", row.key,
                "adsorbate=", row.adsorbate,
                "sum=", row.key_formula,
                "ref=", row.ref,
                "dname=", row.dname,
                "n_formula_matches=", len(matches),
                "formula_matches=", matches,
                "ref_found=", ref_ok,
                "source=", row.source,
            )

    def all_adsorbate_groups(self, row):
        groups = []
        for value in row.adsorbate_dict.values():
            groups.extend(value)
        return groups

    def cross_distance(self, atoms_a, i, atoms_b, j):
        pa = atoms_a.get_positions()[i]
        pb = atoms_b.get_positions()[j]
        vec = pb - pa
        try:
            from ase.geometry import find_mic
            vec, _ = find_mic(vec, cell=atoms_a.get_cell().array, pbc=atoms_a.get_pbc())
        except Exception:
            pass
        return float((vec * vec).sum() ** 0.5)

    def cross_group_distance(self, row_a, group_a, row_b, group_b):
        best = float("inf")
        for i in group_a:
            for j in group_b:
                d = self.cross_distance(row_a.atoms, i, row_b.atoms, j)
                if d < best:
                    best = d
        return best

    def group_distance_same_row(self, row, group_a, group_b):
        best = float("inf")
        for i in group_a:
            for j in group_b:
                d = self.mic_distance(row.atoms, i, j)
                if d < best:
                    best = d
        return best

    def choose_nearest_key_units(self, row):
        """
        Wenn row.key mehrere Einheiten enthält, aber die komplette Summenformel
        nicht als ein Adsorbat vorkommt:

          A+B -> nächstes Paar A und B
          3A  -> nächstes Tripel A,A,A

        Score = Summe aller paarweisen minimalen Abstände.
        Kleiner Score heißt räumlich nächste Kombination.
        """
        unit_formulas = self.expanded_key_unit_formulas(row.key)

        if len(unit_formulas) <= 1:
            return None, "no-exact-formula", {
                "unit_formulas": unit_formulas,
                "combo_scores": [],
                "chosen_score": None,
            }

        candidate_lists = []
        for formula in unit_formulas:
            candidates = row.adsorbate_dict.get(formula, [])
            if not candidates:
                return None, "multiunit-missing-unit-%s" % formula, {
                    "unit_formulas": unit_formulas,
                    "combo_scores": [],
                    "chosen_score": None,
                }
            candidate_lists.append(candidates)

        from itertools import product, combinations

        best_score = None
        best_combo = None
        combo_scores = []
        n_checked = 0

        for combo in product(*candidate_lists):
            # Keine Komponente doppelt verwenden.
            used = set()
            ok = True
            for group in combo:
                group_set = set(group)
                if used & group_set:
                    ok = False
                    break
                used.update(group_set)

            if not ok:
                continue

            score = 0.0
            for a, b in combinations(combo, 2):
                score += self.group_distance_same_row(row, a, b)

            indices = sorted(used)
            combo_scores.append((score, [list(x) for x in combo], indices))
            n_checked += 1

            if best_score is None or score < best_score:
                best_score = score
                best_combo = (combo, indices)

        if best_combo is None:
            return None, "multiunit-no-valid-combination", {
                "unit_formulas": unit_formulas,
                "combo_scores": combo_scores,
                "chosen_score": None,
            }

        combo, indices = best_combo
        return indices, "multiunit-nearest-combination", {
            "unit_formulas": unit_formulas,
            "combo_scores": sorted(combo_scores, key=lambda x: x[0])[:20],
            "chosen_score": best_score,
            "chosen_groups": [list(x) for x in combo],
        }

    def format_combo_scores_for_log(self, combo_scores):
        if not combo_scores:
            return "[]"

        parts = []
        for score, groups, indices in combo_scores:
            parts.append(
                "groups=%s indices=%s score=%.6f" % (groups, indices, score)
            )
        return "[" + "; ".join(parts) + "]"

    def choose_exact_match(self, row, matches):
        """
        Fall: komplette Summenformel existiert als Adsorbat-Komponente.

        Return:
            indices, method, info

        info enthält für das Log u.a. alle Kandidaten und die Scores.
        """
        if len(matches) == 0:
            indices, method, info = self.choose_nearest_key_units(row)
            info["matches"] = matches
            info["scores"] = []
            return indices, method, info

        if len(matches) == 1:
            return matches[0], "exact-single", {
                "matches": matches,
                "scores": [(None, matches[0])],
                "chosen_score": None,
            }

        # Mehrere Adsorbate mit derselben Summenformel:
        # jetzt muss die Referenz entscheiden.
        if row.ref_row is None:
            return None, "same-formula-multiple-no-ref-row", {
                "matches": matches,
                "scores": [],
                "chosen_score": None,
            }

        ref_groups = self.all_adsorbate_groups(row.ref_row)
        if not ref_groups:
            return None, "same-formula-multiple-empty-ref", {
                "matches": matches,
                "scores": [],
                "chosen_score": None,
            }

        scored = []
        for match in matches:
            min_dist = min(
                self.cross_group_distance(row, match, row.ref_row, ref_group)
                for ref_group in ref_groups
            )
            scored.append((min_dist, match))

        scored.sort(key=lambda x: x[0], reverse=True)
        chosen_score, chosen = scored[0]
        return chosen, "same-formula-multiple-farthest-from-ref", {
            "matches": matches,
            "scores": scored,
            "chosen_score": chosen_score,
        }

    def format_scores_for_log(self, scores):
        if not scores:
            return "[]"
        parts = []
        for score, indices in scores:
            if score is None:
                parts.append("%s" % indices)
            else:
                parts.append("%s:min_dist=%.6f" % (indices, score))
        return "[" + "; ".join(parts) + "]"

    def write_details_exact(self, path, log_path=None, debug=False, write_log=False):
        """
        Schreibt details_new.txt.

        Normalmodus:
            - keine Prints
            - kein Logfile
            - nur details_new.txt

        --debug-write:
            - KEEP/SKIP nach stderr

        --write-log:
            - details_new.log mit vollständigem Entscheidungslog

        Wichtig:
            Diese Schleife läuft über die komplette Tabelle:
                for i, row in enumerate(self.rows)
        """
        lines = []
        selected = []
        log_lines = []
        stats = Counter()

        if write_log:
            log_lines.append("# details_new generation log")
            log_lines.append("# output_file: %s" % path)
            log_lines.append("# bond_scale: %s" % self.bond_scale)
            log_lines.append("# ts_bond_scale: %s" % self.ts_bond_scale)
            log_lines.append("# columns in details_new.txt: adsorbate atom_indices ref")
            log_lines.append("#")

        for i, row in enumerate(self.rows):
            matches = row.adsorbate_dict.get(row.key_formula, []) if row.key_formula else []
            indices, method, info = self.choose_exact_match(row, matches)
            stats[method] += 1

            prefix = (
                "ROW %d | method=%s | key=%s | adsorbate=%s | sum=%s | "
                "ref=%s | dname=%s | source=%s | n_formula_matches=%d"
                % (
                    i,
                    method,
                    row.key,
                    row.adsorbate,
                    row.key_formula,
                    row.ref,
                    row.dname,
                    row.source,
                    len(matches),
                )
            )

            if indices is None:
                if write_log:
                    log_lines.append(
                        prefix
                        + " | action=SKIP | ref_found=%s | formula_matches=%s | unit_formulas=%s"
                        % (row.ref_row is not None, matches, info.get("unit_formulas", []))
                    )

                if debug:
                    print(
                        "SKIP_ROW",
                        i,
                        "method=", method,
                        "key=", row.key,
                        "adsorbate=", row.adsorbate,
                        "sum=", row.key_formula,
                        "ref=", row.ref,
                        "dname=", row.dname,
                        "n_formula_matches=", len(matches),
                        "formula_matches=", matches,
                        "ref_found=", row.ref_row is not None,
                        file=sys.stderr,
                    )
                continue

            if not row.ref:
                stats["skip-no-ref"] += 1

                if write_log:
                    log_lines.append(
                        prefix
                        + " | action=SKIP | reason=no-ref | chosen_indices=%s"
                        % indices
                    )

                if debug:
                    print(
                        "SKIP_ROW",
                        i,
                        "method=skip-no-ref",
                        "key=", row.key,
                        "adsorbate=", row.adsorbate,
                        "sum=", row.key_formula,
                        "ref=", row.ref,
                        "dname=", row.dname,
                        "indices=", indices,
                        file=sys.stderr,
                    )
                continue

            line = "%s %s %s" % (row.adsorbate, ",".join(str(x) for x in indices), row.ref)
            lines.append(line)
            selected.append(
                {
                    "row_index": i,
                    "row": row,
                    "indices": list(indices),
                    "method": method,
                    "line": line,
                }
            )

            if write_log:
                if method == "exact-single":
                    log_lines.append(
                        prefix
                        + " | action=KEEP | write_details=yes"
                        + " | reason=unique-formula-match"
                        + " | chosen_indices=%s" % indices
                        + " | line=%s" % line
                    )
                elif method == "same-formula-multiple-farthest-from-ref":
                    log_lines.append(
                        prefix
                        + " | action=KEEP | write_details=yes"
                        + " | reason=multiple-same-formula-farthest-from-ref"
                        + " | ref_dname=%s" % (row.ref_row.dname if row.ref_row is not None else "")
                        + " | candidates=%s" % self.format_scores_for_log(info.get("scores", []))
                        + " | chosen_indices=%s" % indices
                        + " | chosen_min_dist=%s" % (
                            "%.6f" % info.get("chosen_score")
                            if info.get("chosen_score") is not None
                            else "None"
                        )
                        + " | line=%s" % line
                    )
                elif method == "multiunit-nearest-combination":
                    log_lines.append(
                        prefix
                        + " | action=KEEP | write_details=yes"
                        + " | reason=multiunit-nearest-combination"
                        + " | unit_formulas=%s" % info.get("unit_formulas", [])
                        + " | chosen_groups=%s" % info.get("chosen_groups", [])
                        + " | chosen_indices=%s" % indices
                        + " | chosen_score=%s" % (
                            "%.6f" % info.get("chosen_score")
                            if info.get("chosen_score") is not None
                            else "None"
                        )
                        + " | candidates=%s" % self.format_combo_scores_for_log(info.get("combo_scores", []))
                        + " | line=%s" % line
                    )
                else:
                    log_lines.append(
                        prefix
                        + " | action=KEEP | write_details=yes"
                        + " | reason=%s" % method
                        + " | chosen_indices=%s" % indices
                        + " | line=%s" % line
                    )

            if debug:
                print(
                    "KEEP_ROW",
                    i,
                    "method=", method,
                    "key=", row.key,
                    "adsorbate=", row.adsorbate,
                    "sum=", row.key_formula,
                    "ref=", row.ref,
                    "dname=", row.dname,
                    "indices=", indices,
                    "line=", line,
                    file=sys.stderr,
                )

        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")

        if write_log and log_path is not None:
            log_lines.append("#")
            log_lines.append("# stats")
            for key in sorted(stats):
                log_lines.append("STAT %s %d" % (key, stats[key]))
            log_lines.append("STAT written_lines %d" % len(lines))

            with open(log_path, "w") as f:
                for line in log_lines:
                    f.write(line + "\n")

        return stats, len(lines), selected

    def write_selected_adsorbate_traj(self, selected, traj_path):
        """
        Schreibt die ausgewählten Adsorbate als isolierte Moleküle in eine .traj.

        Jede Struktur enthält nur die gewählten Atomindizes aus details_new.txt.
        Metadata liegt in atoms.info:
            adsorbate, ref, source_dname, method, row_index, details_line,
            original_indices
        """
        if traj_path is None:
            return 0

        try:
            from ase.io import write
        except Exception as exc:
            raise RuntimeError(
                "ASE konnte nicht importiert werden; .traj-Ausgabe nicht möglich"
            ) from exc

        molecules = []

        for item in selected:
            row = item["row"]
            indices = list(item["indices"])

            if row.atoms is None or not indices:
                continue

            mol = row.atoms[indices].copy()
            mol.info["adsorbate"] = row.adsorbate
            mol.info["ref"] = row.ref
            mol.info["source_dname"] = row.dname
            mol.info["method"] = item["method"]
            mol.info["row_index"] = item["row_index"]
            mol.info["details_line"] = item["line"]
            mol.info["original_indices"] = ",".join(str(x) for x in indices)

            molecules.append(mol)

        write(str(traj_path), molecules)
        return len(molecules)

def backup_existing_file(path):
    """
    Sicherheitsnetz: bevor details_new.txt überschrieben wird, wird eine
    fortlaufende Backup-Datei geschrieben.

    Beispiel:
        details_new.txt.bak1
        details_new.txt.bak2
    """
    if path is None or not path.exists():
        return None

    n = 1
    while True:
        backup = Path(str(path) + ".bak%d" % n)
        if not backup.exists():
            backup.write_text(path.read_text())
            return backup
        n += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outputs_pickle", type=Path, default=Path("structures.pkl"))
    parser.add_argument("--compounds_pickle", type=Path, default=Path("compounds.pkl"))
    parser.add_argument("--surface-elements", default=None)
    parser.add_argument("--bond-scale", type=float, default=1.35)
    parser.add_argument("--ts-bond-scale", type=float, default=1.75)
    parser.add_argument("-o", "--output", type=Path, default=Path("details_new.txt"))
    parser.add_argument("--no-write-details", action="store_true")
    parser.add_argument("--no-backup", action="store_true")

    # Logfile nur schreiben, wenn diese Option gesetzt ist.
    parser.add_argument("--write-log", action="store_true")
    parser.add_argument("--log", type=Path, default=Path("details_new.log"))

    # Optional: ausgewählte Adsorbat-Fragmente als Moleküle in eine .traj schreiben.
    # Ohne Dateiname wird adsorbates.traj verwendet.
    parser.add_argument(
        "--adsorbate-traj",
        nargs="?",
        const=Path("adsorbates.traj"),
        default=None,
        type=Path,
    )

    parser.add_argument("--print-table", action="store_true")
    parser.add_argument("--test-exact", action="store_true")
    parser.add_argument("--debug-write", action="store_true")

    args = parser.parse_args()

    data = CompoundData.from_files(
        args.outputs_pickle,
        args.compounds_pickle,
        surface_elements=args.surface_elements,
        bond_scale=args.bond_scale,
        ts_bond_scale=args.ts_bond_scale,
    )

    if args.print_table:
        data.print_table()

    if args.test_exact:
        data.print_exact_formula_test()

    # Sicherheitslogik:
    # - normaler Lauf schreibt details_new.txt
    # - --print-table / --test-exact allein dienen zum Debuggen und schreiben NICHT
    #   automatisch details_new.txt
    # - --no-write-details verhindert das Schreiben explizit
    debug_only_mode = (args.print_table or args.test_exact) and not args.adsorbate_traj and not args.write_log
    write_details = not args.no_write_details and not debug_only_mode

    stats = Counter()
    n_lines = 0
    selected = []
    backup = None

    if write_details:
        if not args.no_backup:
            backup = backup_existing_file(args.output)

        stats, n_lines, selected = data.write_details_exact(
            args.output,
            log_path=args.log,
            debug=args.debug_write,
            write_log=args.write_log,
        )
    else:
        # Für reine Debug-/Print-Läufe kein details_new.txt anfassen.
        # Falls nur eine Trajectory gewünscht ist, brauchen wir trotzdem die
        # Auswahl; dafür schreiben wir temporär NICHT nach details_new.txt,
        # sondern bauen die Auswahl über write_details_exact in eine temp-Datei.
        if args.adsorbate_traj is not None:
            tmp_output = Path(str(args.output) + ".tmp_for_traj")
            stats, n_lines, selected = data.write_details_exact(
                tmp_output,
                log_path=args.log,
                debug=args.debug_write,
                write_log=args.write_log,
            )
            try:
                tmp_output.unlink()
            except FileNotFoundError:
                pass

    n_mols = 0
    if args.adsorbate_traj is not None:
        n_mols = data.write_selected_adsorbate_traj(selected, args.adsorbate_traj)

    if args.debug_write:
        print("N_ROWS", len(data.rows), file=sys.stderr)
        print("WRITE_DETAILS", write_details, file=sys.stderr)
        if backup is not None:
            print("BACKUP", backup, file=sys.stderr)
        if write_details:
            print("DETAILS", args.output, file=sys.stderr)
        if args.write_log:
            print("LOG", args.log, file=sys.stderr)
        if args.adsorbate_traj is not None:
            print("ADSORBATE_TRAJ", args.adsorbate_traj, file=sys.stderr)
            print("N_ADSORBATE_MOLS", n_mols, file=sys.stderr)
        print("N_DETAILS", n_lines, file=sys.stderr)
        for key in sorted(stats):
            print("STAT", key, stats[key], file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
