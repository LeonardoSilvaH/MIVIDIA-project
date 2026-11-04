#!/usr/bin/env python3

"""
02_match_pockets.py

Match WT and mutant CoBDock pockets so that docking scores are only
compared within the SAME structural pocket.

Pocket identity
---------------
Two pockets are the same site when they share enough pocket-lining
residues (Jaccard overlap of residue numbers >= MIN_RESIDUE_OVERLAP).
Residues are compared by (chain, residue number) only, so the mutated
residue still counts as shared.

Residue lists come from fpocket_data.json / p2rank_data.json, using the
selected_pocket_program and selected_pocket_id stored for every pocket
in top_3_pockets.json. When the residues of a pocket cannot be found
(or the lookup does not agree with the pocket center), the match falls
back to the aligned center distance (<= MAX_CENTER_DISTANCE).

Outcome for every WT top-3 pocket
---------------------------------
Matched      : the same pocket is among the mutant top-3 (docked) pockets
NotDocked    : the pocket still exists in the mutant, but was not docked
Lost         : no mutant pocket overlaps it
Unresolved   : no docked match, and no residue data to search further

Mutant top-3 pockets without a WT counterpart are reported as NewInMutant.
Only Matched rows are scored by 03_compute_delta_scores.py.

Inputs
------
Analysis/01_align_structures/alignment_matrices.pkl
<results>/pocket_locations/top_3_pockets.json
<results>/pocket_locations/pose_data.json
<results>/docking/fpocket_data.json
<results>/docking/p2rank_data.json

Outputs
-------
Analysis/02_matched_pockets/02_matched_pockets.csv
Analysis/02_matched_pockets/02_distance_matrix.csv
Analysis/02_matched_pockets/02_skipped_entries.csv
"""

from pathlib import Path
from dataclasses import dataclass
import json
import pickle

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


###############################################################################
# USER SETTINGS
###############################################################################

ALIGNMENT_FILE = Path(
    "Analysis/01_align_structures/alignment_matrices.pkl"
)

# ---------------------------------------------------------------------------
# CoBDock results
# ---------------------------------------------------------------------------

POCKET_LOCATIONS_DIR = Path(
    "/hd/leo/MIVIDIA/Data/resultados_mestrado/pocket_locations"
)

POSE_JSON = (
    POCKET_LOCATIONS_DIR /
    "pose_data.json"
)

POCKET_JSON = (
    POCKET_LOCATIONS_DIR /
    "top_3_pockets.json"
)

# Pocket residue lists, written by CoBDock next to pocket_locations/
DOCKING_DIR = (
    POCKET_LOCATIONS_DIR.parent /
    "docking"
)

POCKET_RESIDUE_FILES = {
    "fpocket": DOCKING_DIR / "fpocket_data.json",
    "p2rank": DOCKING_DIR / "p2rank_data.json",
}

# ---------------------------------------------------------------------------
# Same-pocket criteria
# ---------------------------------------------------------------------------

# Minimum Jaccard overlap of pocket-lining residues
MIN_RESIDUE_OVERLAP = 0.5

# Maximum aligned center distance (Å), used only when residues are missing
MAX_CENTER_DISTANCE = 4.0

# A residue lookup is trusted only if its center agrees within this (Å)
CENTER_TOLERANCE = 0.01

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path(
    "Analysis/02_matched_pockets"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

MATCHED_OUTPUT = (
    OUTPUT_DIR /
    "02_matched_pockets.csv"
)

DISTANCE_OUTPUT = (
    OUTPUT_DIR /
    "02_distance_matrix.csv"
)

SKIPPED_OUTPUT = (
    OUTPUT_DIR /
    "02_skipped_entries.csv"
)


###############################################################################
# DATA STRUCTURE
###############################################################################

@dataclass
class Pocket:

    rank: int
    program: str
    pocket_id: str
    pocket_score: float
    center: np.ndarray
    residues: frozenset = None
    docking_score: float = np.nan


###############################################################################
# LOAD FILES
###############################################################################

def load_json(path):

    with open(path) as f:

        return json.load(f)


print("Loading alignment matrices...")

with open(ALIGNMENT_FILE, "rb") as f:

    alignments = pickle.load(f)

print("Loading top_3_pockets.json...")

pocket_json = load_json(POCKET_JSON)

print("Loading pose_data.json...")

pose_json = load_json(POSE_JSON)

pocket_residue_data = {}

for program, path in POCKET_RESIDUE_FILES.items():

    if path.exists():

        print(f"Loading {path.name}...")

        pocket_residue_data[program] = load_json(path)

    else:

        print(
            f"WARNING: {path} not found. {program} pockets will be "
            f"matched by center distance only."
        )

        pocket_residue_data[program] = {}


###############################################################################
# HELPER FUNCTIONS
###############################################################################

def get_wildtype_key(entry):
    """
    Find the WT entry.

    Examples:

        WILD-TYE
        KRAS-WILD-TYE
        KRAS.WILD-TYE
    """

    for key in entry.keys():

        if "WILD" in key.upper():

            return key

    raise ValueError(
        "Wild-type not found."
    )


def alignment_mutation_name(protein, mutation):
    """
    Convert a CoBDock mutation key to the name used by
    01_align_structures.py ("PROTEIN-MUTATION").

        ACVR1.G328E  ->  ACVR1-G328E
        G12D         ->  KRAS-G12D
        KRAS-G12D    ->  KRAS-G12D
    """

    for separator in (".", "-"):

        prefix = f"{protein}{separator}"

        if mutation.startswith(prefix):

            return f"{protein}-" + mutation[len(prefix):]

    return f"{protein}-{mutation}"


def split_ligand_name(name):
    """
    Split:

        MYC_Prednisone

    into:

        MYC
        Prednisone
    """

    protein, drug = name.split(
        "_",
        1
    )

    return protein, drug


def pocket_center(pocket):

    return np.array(
        [
            pocket["center_x"],
            pocket["center_y"],
            pocket["center_z"]
        ],
        dtype=float
    )


def residue_set(residue_ids):
    """
    [["A", "GLY", 12], ...]  ->  {("A", 12), ...}

    The residue name is ignored, so the mutated residue is still
    recognised as the same position in WT and mutant.
    """

    residues = set()

    for residue in residue_ids or []:

        try:

            residues.add(
                (
                    str(residue[0]),
                    int(residue[-1])
                )
            )

        except (
            IndexError,
            TypeError,
            ValueError
        ):

            continue

    return frozenset(residues) or None


def lookup_residues(program, mutation, structure, pocket_id, center):
    """
    Pocket-lining residues of one CoBDock pocket, or None when they
    cannot be found or belong to a different pocket (center mismatch).
    """

    pocket = (
        pocket_residue_data
        .get(program, {})
        .get(mutation, {})
        .get(structure, {})
        .get(str(pocket_id))
    )

    if pocket is None:

        return None

    try:

        lookup_center = pocket_center(pocket)

    except (
        KeyError,
        TypeError,
        ValueError
    ):

        return None

    if np.abs(lookup_center - center).max() > CENTER_TOLERANCE:

        return None

    return residue_set(
        pocket.get("residue_ids")
    )


def best_docking_scores(structure_dict):
    """
    Extract the best PLANTS score for each pocket rank.
    """

    scores = {}

    for pocket_number, poses in structure_dict.items():

        best = None

        for pose in poses.values():

            try:

                score = float(
                    pose.get("total_score")
                )

            except (
                TypeError,
                ValueError
            ):

                continue

            if (
                best is None
                or score < best
            ):

                best = score

        scores[int(pocket_number)] = (
            np.nan if best is None else best
        )

    return scores


def build_docked_pockets(
    pocket_list,
    docking_scores,
    mutation,
    structure
):
    """
    Build Pocket objects for the top-3 (docked) pockets.
    """

    pockets = []

    for pocket in pocket_list:

        rank = int(pocket["pocket_rank"])

        program = str(
            pocket.get("selected_pocket_program", "fpocket")
        )

        pocket_id = str(pocket["selected_pocket_id"])

        center = pocket_center(pocket)

        pockets.append(
            Pocket(
                rank=rank,
                program=program,
                pocket_id=pocket_id,
                pocket_score=float(pocket["pocket_score"]),
                center=center,
                residues=lookup_residues(
                    program,
                    mutation,
                    structure,
                    pocket_id,
                    center
                ),
                docking_score=docking_scores.get(rank, np.nan),
            )
        )

    pockets.sort(
        key=lambda x: x.rank
    )

    return pockets


def all_detected_pockets(mutation, structure):
    """
    Every pocket detected by fpocket and p2rank in one structure,
    including those that were not docked.
    """

    pockets = []

    for program, data in pocket_residue_data.items():

        structure_pockets = (
            data
            .get(mutation, {})
            .get(structure, {})
        )

        for pocket_id, pocket in structure_pockets.items():

            residues = residue_set(
                pocket.get("residue_ids")
            )

            if residues is None:

                continue

            try:

                center = pocket_center(pocket)

            except (
                KeyError,
                TypeError,
                ValueError
            ):

                continue

            pockets.append(
                Pocket(
                    rank=None,
                    program=program,
                    pocket_id=str(pocket_id),
                    pocket_score=float(
                        pocket.get("pocket_score", np.nan)
                    ),
                    center=center,
                    residues=residues,
                )
            )

    return pockets


def residue_overlap(a, b):
    """
    Jaccard overlap of two residue sets, or None if either is missing.
    """

    if a is None or b is None:

        return None

    return len(a & b) / len(a | b)


def transform_center(center, alignment):
    """
    Move a mutant pocket center into the WT reference frame.
    """

    if alignment is None:

        return None

    return (
        center @ alignment["rotation"]
    ) + alignment["translation"]


def compare_pockets(wt, mut, alignment):
    """
    Decide whether a WT and a mutant pocket are the same site.

    Returns
    -------
    dict with ResidueOverlap, Distance, MatchMethod, SamePocket
    """

    overlap = residue_overlap(
        wt.residues,
        mut.residues
    )

    mut_center = transform_center(
        mut.center,
        alignment
    )

    distance = (
        np.nan if mut_center is None
        else float(np.linalg.norm(wt.center - mut_center))
    )

    if overlap is not None:

        method = "Residues"

        same = overlap >= MIN_RESIDUE_OVERLAP

    elif not np.isnan(distance):

        method = "Distance"

        same = distance <= MAX_CENTER_DISTANCE

    else:

        method = "None"

        same = False

    return {
        "ResidueOverlap": np.nan if overlap is None else overlap,
        "Distance": distance,
        "MatchMethod": method,
        "SamePocket": same,
        "MutCenterAligned": mut_center,
    }


def match_cost(comparison):
    """
    Lower is better. Residue overlap is preferred; distance only
    decides when residues are unavailable.
    """

    if not comparison["SamePocket"]:

        return 1e6

    if comparison["MatchMethod"] == "Residues":

        return 1.0 - comparison["ResidueOverlap"]

    return 1.0 + comparison["Distance"] / MAX_CENTER_DISTANCE


###############################################################################
# MAIN ANALYSIS
###############################################################################

matched_rows = []

distance_rows = []

skipped_rows = []


def skip_entry(protein, drug, mutation, reason):
    """
    Record an entry that could not be processed instead of crashing.
    """

    print(
        f"    Skipping {protein} - {drug} - {mutation}: {reason}"
    )

    skipped_rows.append(
        {
            "Protein": protein,
            "Drug": drug,
            "Mutation": mutation,
            "Reason": reason,
        }
    )


def pocket_fields(prefix, pocket, center=None):

    if pocket is None:

        return {
            f"{prefix}_PocketRank": np.nan,
            f"{prefix}_PocketProgram": None,
            f"{prefix}_PocketID": None,
            f"{prefix}_PocketScore": np.nan,
            f"{prefix}_DockingScore": np.nan,
            f"{prefix}_NResidues": np.nan,
            f"{prefix}_X": np.nan,
            f"{prefix}_Y": np.nan,
            f"{prefix}_Z": np.nan,
        }

    if center is None:

        center = np.full(3, np.nan)

    return {
        f"{prefix}_PocketRank": pocket.rank,
        f"{prefix}_PocketProgram": pocket.program,
        f"{prefix}_PocketID": pocket.pocket_id,
        f"{prefix}_PocketScore": pocket.pocket_score,
        f"{prefix}_DockingScore": pocket.docking_score,
        f"{prefix}_NResidues": (
            np.nan if pocket.residues is None
            else len(pocket.residues)
        ),
        f"{prefix}_X": center[0],
        f"{prefix}_Y": center[1],
        f"{prefix}_Z": center[2],
    }


def match_row(protein, drug, mutation, status, wt, mut, comparison):

    comparison = comparison or {}

    return {
        "Protein": protein,
        "Drug": drug,
        "Mutation": mutation,
        "MatchStatus": status,
        "MatchMethod": comparison.get("MatchMethod"),
        "ResidueOverlap": comparison.get("ResidueOverlap", np.nan),
        "Distance": comparison.get("Distance", np.nan),
        **pocket_fields(
            "WT",
            wt,
            None if wt is None else wt.center
        ),
        # Mutant coordinates are given in the WT reference frame
        **pocket_fields(
            "Mut",
            mut,
            comparison.get("MutCenterAligned")
        ),
    }


for ligand_name in sorted(pocket_json):

    protein, drug = split_ligand_name(
        ligand_name
    )

    print(
        f"\nProcessing {protein} - {drug}"
    )

    pocket_entry = pocket_json[ligand_name]

    if ligand_name not in pose_json:

        skip_entry(
            protein,
            drug,
            "ALL",
            "Missing in pose_data.json"
        )

        continue

    pose_entry = pose_json[ligand_name]

    # -------------------------------------------------------------------------
    # WT
    # -------------------------------------------------------------------------

    try:

        wt_key = get_wildtype_key(
            pocket_entry
        )

    except ValueError:

        skip_entry(
            protein,
            drug,
            "ALL",
            "Wild-type missing in top_3_pockets.json"
        )

        continue

    wt_structure = next(
        iter(pocket_entry[wt_key])
    )

    if wt_structure not in pose_entry.get(wt_key, {}):

        skip_entry(
            protein,
            drug,
            "ALL",
            f"Wild-type ({wt_key}/{wt_structure}) "
            f"missing in pose_data.json"
        )

        continue

    wt_pockets = build_docked_pockets(
        pocket_entry[wt_key][wt_structure],
        best_docking_scores(
            pose_entry[wt_key][wt_structure]
        ),
        wt_key,
        wt_structure
    )

    # -------------------------------------------------------------------------
    # Mutants
    # -------------------------------------------------------------------------

    for mutation in pocket_entry:

        if mutation == wt_key:

            continue

        mut_structure = next(
            iter(pocket_entry[mutation])
        )

        if mut_structure not in pose_entry.get(mutation, {}):

            skip_entry(
                protein,
                drug,
                mutation,
                f"Mutant ({mutation}/{mut_structure}) "
                f"missing in pose_data.json"
            )

            continue

        # Residue overlap does not need the alignment; it is only
        # used for center distances.
        alignment = alignments.get(
            (
                protein,
                alignment_mutation_name(protein, mutation)
            )
        )

        mut_pockets = build_docked_pockets(
            pocket_entry[mutation][mut_structure],
            best_docking_scores(
                pose_entry[mutation][mut_structure]
            ),
            mutation,
            mut_structure
        )

        if alignment is None and any(
            p.residues is None
            for p in wt_pockets + mut_pockets
        ):

            skip_entry(
                protein,
                drug,
                mutation,
                "Alignment missing in alignment_matrices.pkl and "
                "pocket residues unavailable"
            )

            continue

        # ---------------------------------------------------------------------
        # Compare every WT pocket with every docked mutant pocket
        # ---------------------------------------------------------------------

        comparisons = [
            [
                compare_pockets(wt, mut, alignment)
                for mut in mut_pockets
            ]
            for wt in wt_pockets
        ]

        for i, wt in enumerate(wt_pockets):

            for j, mut in enumerate(mut_pockets):

                c = comparisons[i][j]

                distance_rows.append(
                    {
                        "Protein": protein,
                        "Drug": drug,
                        "Mutation": mutation,
                        "WT_PocketRank": wt.rank,
                        "WT_PocketID": wt.pocket_id,
                        "Mut_PocketRank": mut.rank,
                        "Mut_PocketID": mut.pocket_id,
                        "ResidueOverlap": c["ResidueOverlap"],
                        "Distance": c["Distance"],
                        "MatchMethod": c["MatchMethod"],
                        "SamePocket": c["SamePocket"],
                    }
                )

        # ---------------------------------------------------------------------
        # One-to-one assignment, keeping only same-pocket pairs
        # ---------------------------------------------------------------------

        cost = np.array(
            [
                [match_cost(c) for c in row]
                for row in comparisons
            ]
        ).reshape(len(wt_pockets), len(mut_pockets))

        matches = {}

        if cost.size:

            for i, j in zip(*linear_sum_assignment(cost)):

                if comparisons[i][j]["SamePocket"]:

                    matches[i] = j

        # ---------------------------------------------------------------------
        # Outcome of every WT pocket
        # ---------------------------------------------------------------------

        candidates = None

        for i, wt in enumerate(wt_pockets):

            if i in matches:

                j = matches[i]

                matched_rows.append(
                    match_row(
                        protein,
                        drug,
                        mutation,
                        "Matched",
                        wt,
                        mut_pockets[j],
                        comparisons[i][j]
                    )
                )

                continue

            if wt.residues is None:

                matched_rows.append(
                    match_row(
                        protein,
                        drug,
                        mutation,
                        "Unresolved",
                        wt,
                        None,
                        None
                    )
                )

                continue

            # Search every pocket detected in the mutant, docked or not
            if candidates is None:

                candidates = all_detected_pockets(
                    mutation,
                    mut_structure
                )

            best, best_comparison = None, None

            for candidate in candidates:

                c = compare_pockets(wt, candidate, alignment)

                if c["SamePocket"] and (
                    best_comparison is None
                    or c["ResidueOverlap"] > best_comparison["ResidueOverlap"]
                ):

                    best, best_comparison = candidate, c

            matched_rows.append(
                match_row(
                    protein,
                    drug,
                    mutation,
                    "Lost" if best is None else "NotDocked",
                    wt,
                    best,
                    best_comparison
                )
            )

        # ---------------------------------------------------------------------
        # Docked mutant pockets with no WT counterpart
        # ---------------------------------------------------------------------

        used = set(matches.values())

        for j, mut in enumerate(mut_pockets):

            if j in used:

                continue

            matched_rows.append(
                match_row(
                    protein,
                    drug,
                    mutation,
                    "NewInMutant",
                    None,
                    mut,
                    {
                        "MutCenterAligned": transform_center(
                            mut.center,
                            alignment
                        )
                    }
                )
            )


###############################################################################
# CREATE DATAFRAMES
###############################################################################

matched_df = pd.DataFrame(
    matched_rows
)

distance_df = pd.DataFrame(
    distance_rows
)

skipped_df = pd.DataFrame(
    skipped_rows,
    columns=["Protein", "Drug", "Mutation", "Reason"]
)

if not matched_df.empty:

    for column in ("WT_PocketRank", "Mut_PocketRank"):

        matched_df[column] = matched_df[column].astype("Int64")

    matched_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank",
            "Mut_PocketRank"
        ],
        inplace=True
    )

if not distance_df.empty:

    distance_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank",
            "Mut_PocketRank"
        ],
        inplace=True
    )


###############################################################################
# SAVE
###############################################################################

matched_df.to_csv(
    MATCHED_OUTPUT,
    index=False
)

distance_df.to_csv(
    DISTANCE_OUTPUT,
    index=False
)

skipped_df.to_csv(
    SKIPPED_OUTPUT,
    index=False
)


###############################################################################
# SUMMARY
###############################################################################

print()
print("=" * 70)
print("Finished.")
print("=" * 70)

print(f"Pocket rows      : {len(matched_df)}")
print(f"Distance entries : {len(distance_df)}")
print(f"Skipped entries  : {len(skipped_df)}")

print()
print("Matched output:")
print(f"    {MATCHED_OUTPUT}")
print("Distance output:")
print(f"    {DISTANCE_OUTPUT}")
print("Skipped entries log:")
print(f"    {SKIPPED_OUTPUT}")

print("=" * 70)

if matched_df.empty:

    print(
        "\nNo pockets were processed, so no statistics were computed."
    )

    print(
        f"Check {SKIPPED_OUTPUT} to see why each entry was skipped."
    )

else:

    print("\nWT pocket outcome")
    print("-----------------")

    print(
        matched_df["MatchStatus"]
        .value_counts()
        .to_string()
    )

    matched = matched_df[
        matched_df["MatchStatus"] == "Matched"
    ]

    print("\nMatch method (Matched only)")
    print("---------------------------")

    print(
        matched["MatchMethod"]
        .value_counts()
        .to_string()
    )

    overlap = matched["ResidueOverlap"].dropna()

    dist = matched["Distance"].dropna()

    if len(overlap):

        print("\nResidue overlap (Matched)")
        print("-------------------------")
        print(f"Mean:   {overlap.mean():.2f}")
        print(f"Median: {overlap.median():.2f}")
        print(f"Min:    {overlap.min():.2f}")

    if len(dist):

        print("\nCenter distance (Matched)")
        print("-------------------------")
        print(f"Mean:   {dist.mean():.2f} Å")
        print(f"Median: {dist.median():.2f} Å")
        print(f"Max:    {dist.max():.2f} Å")
