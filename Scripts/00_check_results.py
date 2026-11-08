#!/usr/bin/env python3

"""
00_check_results.py

Vet a CoBDock results folder before running the pipeline on it.

Run this after every new CoBDock upload. It is read-only and takes a
few seconds, and it catches the three failures that have each cost a
full pipeline run:

  * the results folder has no docking/ directory, so every pocket is
    matched by centre distance alone
  * the mutation keys use a different separator than the alignments
    ("ACVR1.G328E" vs "ACVR1-G328E"), so no mutant is ever found
  * wild-type entries are shared between proteins, so a protein is
    compared against another protein's pockets

Usage
-----
    python Scripts/00_check_results.py
    python Scripts/00_check_results.py /path/to/results
    MIVIDIA_RESULTS=/path/to/results python Scripts/00_check_results.py

Exit status is 1 if any check fails, so it can gate a batch run.
"""

from pathlib import Path
from collections import defaultdict
import json
import os
import sys


###############################################################################
# SETTINGS
###############################################################################

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_RESULTS = PROJECT_ROOT / "Data" / "projeto_de_mestrado"

# Must agree with 02_match_pockets.py
CENTER_TOLERANCE = 0.01

MIN_RESIDUE_OVERLAP = 0.5

# Residue count below which the Jaccard overlap becomes coarse
SMALL_POCKET = 8

# Approximate canonical lengths, used only to flag impossible residue
# numbers in the wild-type pockets.
PROTEIN_LENGTH = {
    "KRAS": 189,
    "HRAS": 189,
    "NRAS": 189,
    "RRAS2": 204,
    "TP53": 393,
    "IDH1": 414,
    "IDH2": 452,
    "CDK4": 303,
    "GNA11": 359,
    "MYC": 454,
    "AKT1": 480,
    "MAP2K1": 393,
}


###############################################################################
# HELPERS
###############################################################################

def banner(number, title):

    print()
    print("=" * 74)
    print(f"{number}. {title}")
    print("=" * 74)


def find_file(root, name):
    """
    Locate a results file anywhere under root, so a changed folder
    layout does not hide it.
    """

    hits = sorted(root.rglob(name))

    return hits[0] if hits else None


def load(path):

    with open(path) as f:

        return json.load(f)


def residue_key(residue):
    """
    ["A", "GLY", 12] -> ("A", 12). The residue name is dropped so a
    mutated position still counts as the same residue.
    """

    return (str(residue[0]), int(residue[-1]))


###############################################################################
# CHECKS
###############################################################################

def check_layout(root):

    banner(1, "LAYOUT")

    files = {
        name: find_file(root, name)
        for name in (
            "top_3_pockets.json",
            "pose_data.json",
            "fpocket_data.json",
            "p2rank_data.json",
        )
    }

    for name, path in files.items():

        if path is None:

            print(f"  MISSING  {name}")

        else:

            size = path.stat().st_size / 1e6

            print(f"  OK       {name:22s} {size:8.1f} MB  "
                  f"{path.relative_to(root)}")

    return files


def check_keys(pockets):

    banner(2, "KEY NAMING")

    separators = defaultdict(int)

    wt_keys = set()

    for muts in pockets.values():

        for mutation in muts:

            if "WILD" in mutation.upper():

                wt_keys.add(mutation)

            elif "." in mutation:

                separators["dot"] += 1

            elif "-" in mutation:

                separators["hyphen"] += 1

            else:

                separators["bare"] += 1

    pair = next(iter(pockets))

    print(f"  pairs                : {len(pockets)}")
    print(f"  example pair         : {pair!r}")
    print(f"  example mutations    : {list(pockets[pair])[:4]}")
    print(f"  mutation separators  : {dict(separators)}")
    print(f"  wild-type spellings  : {sorted(wt_keys)[:4]}"
          f"{' ...' if len(wt_keys) > 4 else ''}")

    if len(separators) > 1:

        print("\n  WARNING: mixed separators. 02_match_pockets.py normalises")
        print("  both via alignment_mutation_name(), but check the alignments.")

    return True


def check_residues(pockets, residue_data):
    """
    Walk exactly the pockets 02_match_pockets.py selects and confirm
    each one resolves to a usable residue list.
    """

    banner(3, "RESIDUE LOOKUP")

    if not any(residue_data.values()):

        print("  NO residue data. 02_match_pockets.py will refuse to run")
        print("  (override with MIVIDIA_ALLOW_DISTANCE_ONLY=1).")

        return False

    stats = defaultdict(
        lambda: {"n": 0, "missing": 0, "empty": 0, "rejected": 0, "sizes": []}
    )

    worst_gap = 0.0

    for muts in pockets.values():

        for mutation, structures in muts.items():

            for structure, pocket_list in structures.items():

                for pocket in pocket_list:

                    program = str(
                        pocket.get("selected_pocket_program", "fpocket")
                    )

                    pocket_id = str(pocket["selected_pocket_id"])

                    entry = (
                        residue_data
                        .get(program, {})
                        .get(mutation, {})
                        .get(structure, {})
                        .get(pocket_id)
                    )

                    stat = stats[program]

                    stat["n"] += 1

                    if entry is None:

                        stat["missing"] += 1

                        continue

                    try:

                        gap = max(
                            abs(
                                float(entry[f"center_{axis}"])
                                - float(pocket[f"center_{axis}"])
                            )
                            for axis in "xyz"
                        )

                        worst_gap = max(worst_gap, gap)

                        if gap > CENTER_TOLERANCE:

                            stat["rejected"] += 1

                    except (KeyError, TypeError, ValueError):

                        pass

                    residues = entry.get("residue_ids")

                    if residues:

                        stat["sizes"].append(len(residues))

                    else:

                        stat["empty"] += 1

    ok = True

    for program, stat in sorted(stats.items()):

        sizes = sorted(stat["sizes"])

        print(f"\n  {program}: {stat['n']} selected pockets")
        print(f"    not found in residue file : {stat['missing']}")
        print(f"    empty residue list        : {stat['empty']}")
        print(f"    centre-rejected (>{CENTER_TOLERANCE} A) : "
              f"{stat['rejected']}")

        if stat["missing"] or stat["empty"] or stat["rejected"]:

            ok = False

        if sizes:

            small = sum(1 for n in sizes if n <= SMALL_POCKET)

            print(f"    residues per pocket       : min {sizes[0]}  "
                  f"median {sizes[len(sizes) // 2]}  max {sizes[-1]}")
            print(f"    <= {SMALL_POCKET} residues          : {small} "
                  f"({100 * small / len(sizes):.0f}%)")

            if small:

                step = 1.0 / SMALL_POCKET

                print(f"      on these, one residue moves the Jaccard overlap")
                print(f"      by ~{step:.2f}, so MIN_RESIDUE_OVERLAP="
                      f"{MIN_RESIDUE_OVERLAP} is coarse")

    print(f"\n  worst centre disagreement: {worst_gap:.4f} A")

    return ok


def check_wildtype(pockets, residue_data):
    """
    Two independent corruption tests:

      * do different proteins share identical wild-type pockets?
      * do wild-type pockets cite residues beyond the protein's length?
    """

    banner(4, "WILD-TYPE INTEGRITY")

    fingerprints = defaultdict(set)

    max_residue = {}

    for pair, muts in pockets.items():

        protein = pair.split("_", 1)[0]

        for mutation, structures in muts.items():

            if "WILD" not in mutation.upper():

                continue

            for structure, pocket_list in structures.items():

                fingerprints[protein].add(
                    tuple(sorted(
                        (
                            round(float(p["center_x"]), 2),
                            round(float(p["center_y"]), 2),
                            round(float(p["center_z"]), 2),
                        )
                        for p in pocket_list
                    ))
                )

                for data in residue_data.values():

                    entries = (
                        data
                        .get(mutation, {})
                        .get(structure, {})
                    )

                    for entry in entries.values():

                        for residue in (entry.get("residue_ids") or []):

                            try:

                                number = residue_key(residue)[1]

                            except (TypeError, ValueError, IndexError):

                                continue

                            max_residue[protein] = max(
                                max_residue.get(protein, 0),
                                number
                            )

    shared = defaultdict(set)

    for protein, prints in fingerprints.items():

        for fingerprint in prints:

            shared[fingerprint].add(protein)

    collisions = [
        proteins for proteins in shared.values() if len(proteins) > 1
    ]

    print(f"  proteins with wild-type data : {len(fingerprints)}")

    ok = True

    if collisions:

        ok = False

        print(f"\n  *** {len(collisions)} wild-type pockets shared between "
              f"proteins ***")

        for proteins in collisions[:5]:

            print(f"    {sorted(proteins)}")

        print("\n  This is the corrupt-data signature. Do not use this folder.")

    else:

        print("  OK - every protein has distinct wild-type pockets")

    impossible = [
        (protein, number, PROTEIN_LENGTH[protein])
        for protein, number in max_residue.items()
        if protein in PROTEIN_LENGTH and number > PROTEIN_LENGTH[protein]
    ]

    if impossible:

        ok = False

        print("\n  *** residue numbers beyond the protein length ***")

        for protein, number, length in impossible:

            print(f"    {protein}: residue {number} but length ~{length}")

    elif max_residue:

        print("  OK - no residue number exceeds its protein's length")

    return ok


###############################################################################
# MAIN
###############################################################################

def main():

    root = Path(
        sys.argv[1] if len(sys.argv) > 1
        else os.environ.get("MIVIDIA_RESULTS", DEFAULT_RESULTS)
    ).expanduser().resolve()

    print(f"Checking: {root}")

    if not root.is_dir():

        print(f"\nERROR: not a directory: {root}")

        return 1

    files = check_layout(root)

    if files["top_3_pockets.json"] is None:

        print("\nERROR: top_3_pockets.json not found - cannot continue.")

        return 1

    pockets = load(files["top_3_pockets.json"])

    residue_data = {
        program: (
            load(files[f"{program}_data.json"])
            if files[f"{program}_data.json"] else {}
        )
        for program in ("fpocket", "p2rank")
    }

    results = [
        check_keys(pockets),
        check_residues(pockets, residue_data),
        check_wildtype(pockets, residue_data),
    ]

    banner(5, "VERDICT")

    if all(results):

        print("  All checks passed. Safe to run the pipeline on this folder:")
        print(f"\n    export MIVIDIA_RESULTS={root}")

        return 0

    print("  One or more checks FAILED. See the sections above.")

    return 1


if __name__ == "__main__":

    raise SystemExit(main())
