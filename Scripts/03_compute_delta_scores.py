#!/usr/bin/env python3

import json
import pandas as pd
from tqdm import tqdm
from pathlib import Path


# ==========================================================
# INPUT / OUTPUT
# ==========================================================

MATCHED_POCKETS_FILE = (
    "Analysis/02_matched_pockets/02_matched_pockets.csv"
)

# ----------------------------------------------------------
# CoBDock input
# ----------------------------------------------------------

POCKET_LOCATIONS_DIR = Path(
    "/hd/leo/MIVIDIA/Data/resultados_mestrado/pocket_locations"
)

POSE_DATA_FILE = (
    POCKET_LOCATIONS_DIR /
    "pose_data.json"
)

# ----------------------------------------------------------
# Outputs
# ----------------------------------------------------------

OUTPUT_DIR = Path(
    "Analysis/03_compute_delta_scores"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_FULL = (
    OUTPUT_DIR /
    "matched_pockets_scored.csv"
)

OUTPUT_SUMMARY = (
    OUTPUT_DIR /
    "docking_summary.csv"
)


# ==========================================================
# HELPERS
# ==========================================================

def load_json(filename):

    with open(
        filename,
        "r"
    ) as f:

        return json.load(f)


# ==========================================================
# STANDARD PLANTS FUNCTIONS
# ==========================================================

def get_chain_name(mutation_dict):
    """
    Returns the only chain available.
    """

    return next(
        iter(
            mutation_dict.keys()
        )
    )


def get_best_pose(pocket_dict):
    """
    Find the pose with the lowest PLANTS total_score.

    Returns
    -------
    best_pose : str
    best_score : float
    """

    best_pose = None

    best_score = float(
        "inf"
    )

    for pose_id, pose_data in pocket_dict.items():

        # Skip poses without a numeric score (failed or empty docking)
        try:
            score = float(
                pose_data["total_score"]
            )
        except (
            KeyError,
            TypeError,
            ValueError
        ):
            continue

        if score < best_score:

            best_score = score

            best_pose = pose_id

    if best_pose is None:

        raise ValueError(
            "No pose with a valid PLANTS total_score in this pocket "
            f"({len(pocket_dict)} poses found)"
        )

    return (
        best_pose,
        best_score
    )


def get_pocket(
    pose_data,
    protein,
    drug,
    mutation,
    pocket_rank,
):
    """
    Returns the dictionary corresponding
    to one pocket.
    """

    complex_name = (
        f"{protein}_{drug}"
    )

    mutation_dict = (
        pose_data[
            complex_name
        ][
            mutation
        ]
    )

    chain = get_chain_name(
        mutation_dict
    )

    return (
        mutation_dict[
            chain
        ][
            str(
                int(
                    pocket_rank
                )
            )
        ]
    )


# ==========================================================
# PROCESS ONE MATCH — PLANTS
# ==========================================================

def process_match_plants(
    row,
    pose_data
):

    protein = row["Protein"]

    drug = row["Drug"]

    mutation = row["Mutation"]

    # WT key is "KRAS-WILD-TYE" in older runs and
    # "KRAS.WILD-TYE" in newer ones, so look it up.
    wt_mutation = next(
        key
        for key in pose_data[f"{protein}_{drug}"]
        if "WILD" in key.upper()
    )

    wt_rank = row[
        "WT_PocketRank"
    ]

    mut_rank = row[
        "Mut_PocketRank"
    ]

    # ------------------------------------------------------
    # WT pocket
    # ------------------------------------------------------

    wt_pocket = get_pocket(
        pose_data,
        protein,
        drug,
        wt_mutation,
        wt_rank,
    )

    # ------------------------------------------------------
    # Mutant pocket
    # ------------------------------------------------------

    mut_pocket = get_pocket(
        pose_data,
        protein,
        drug,
        mutation,
        mut_rank,
    )

    # ------------------------------------------------------
    # Best poses
    # ------------------------------------------------------

    wt_pose, wt_score = (
        get_best_pose(
            wt_pocket
        )
    )

    mut_pose, mut_score = (
        get_best_pose(
            mut_pocket
        )
    )

    # ------------------------------------------------------
    # Store scores
    # ------------------------------------------------------

    row["WT_BestPose"] = int(
        wt_pose
    )

    row["WT_BestDockingScore"] = (
        wt_score
    )

    row["Mut_BestPose"] = int(
        mut_pose
    )

    row["Mut_BestDockingScore"] = (
        mut_score
    )

    # ------------------------------------------------------
    # Delta
    # ------------------------------------------------------

    delta = (
        mut_score -
        wt_score
    )

    row["DeltaDockingScore"] = delta

    row["AbsDeltaDockingScore"] = abs(
        delta
    )

    # ------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------

    if delta < 0:

        row["AffinityChange"] = (
            "Stronger"
        )

    elif delta > 0:

        row["AffinityChange"] = (
            "Weaker"
        )

    else:

        row["AffinityChange"] = (
            "Unchanged"
        )

    row["DockingScoreType"] = (
        "PLANTS_PocketSpecific"
    )

    return row


# ==========================================================
# MAIN
# ==========================================================

def main():

    print(
        "Loading data..."
    )

    # ------------------------------------------------------
    # Matched pockets
    # ------------------------------------------------------

    matched_df = pd.read_csv(
        MATCHED_POCKETS_FILE,
        dtype={
            "WT_PocketID": str,
            "Mut_PocketID": str,
        }
    )

    # ------------------------------------------------------
    # Only compare docking scores within the same pocket
    # ------------------------------------------------------

    status_counts = (
        matched_df["MatchStatus"]
        .value_counts()
    )

    print(
        "Pocket outcomes from 02_match_pockets.py:"
    )

    for status, count in status_counts.items():

        print(
            f"  {status:<12} {count}"
        )

    matched_df = matched_df[
        matched_df["MatchStatus"] == "Matched"
    ].copy()

    for column in ("WT_PocketRank", "Mut_PocketRank"):

        matched_df[column] = matched_df[column].astype(int)

    print(
        f"Scoring {len(matched_df)} same-pocket pairs."
    )

    # ------------------------------------------------------
    # Load PLANTS poses
    # ------------------------------------------------------

    print(
        "Loading pose_data.json..."
    )

    docking_data = load_json(
        POSE_DATA_FILE
    )

    # ------------------------------------------------------
    # Process
    # ------------------------------------------------------

    processed_rows = []

    errors = []

    print(
        f"\nProcessing "
        f"{len(matched_df)} matched pockets...\n"
    )

    for index, row in tqdm(
        matched_df.iterrows(),
        total=len(matched_df)
    ):

        try:

            processed_row = (
                process_match_plants(
                    row.copy(),
                    docking_data
                )
            )

            processed_rows.append(
                processed_row
            )

        except Exception as e:

            error_message = (
                f"{row['Protein']} | "
                f"{row['Drug']} | "
                f"{row['Mutation']} | "
                f"{e}"
            )

            errors.append(
                error_message
            )

            print(
                f"\nERROR: {error_message}"
            )

    # ------------------------------------------------------
    # DataFrame
    # ------------------------------------------------------

    scored_df = pd.DataFrame(
        processed_rows
    )

    if scored_df.empty:

        raise RuntimeError(
            "No rows were successfully processed."
        )

    # ------------------------------------------------------
    # Sort
    # ------------------------------------------------------

    scored_df = scored_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank"
        ]
    )

    # ------------------------------------------------------
    # Save full table
    # ------------------------------------------------------

    print(
        "\nSaving complete table..."
    )

    scored_df.to_csv(
        OUTPUT_FULL,
        index=False
    )

    # ======================================================
    # SUMMARY
    # ======================================================

    print(
        "Generating summary table..."
    )

    # Best-scoring same pocket per protein / drug / mutation
    summary_df = (
        scored_df
        .sort_values(
            "Mut_BestDockingScore"
        )
        .groupby(
            [
                "Protein",
                "Drug",
                "Mutation"
            ],
            as_index=False
        )
        .first()
    )

    summary_df.to_csv(
        OUTPUT_SUMMARY,
        index=False
    )

    # ======================================================
    # FINAL REPORT
    # ======================================================

    print()
    print(
        "=" * 60
    )

    print(
        "DONE"
    )

    print(
        "=" * 60
    )

    print(
        f"Same-pocket rows   : {len(matched_df)}"
    )

    print(
        f"Processed rows    : {len(scored_df)}"
    )

    print(
        f"Errors            : {len(errors)}"
    )

    print()

    print(
        f"Full output:"
    )

    print(
        f"  {OUTPUT_FULL}"
    )

    print()

    print(
        f"Summary output:"
    )

    print(
        f"  {OUTPUT_SUMMARY}"
    )

    print(
        "=" * 60
    )

    # ------------------------------------------------------
    # Print errors
    # ------------------------------------------------------

    if errors:

        print()
        print(
            "Errors encountered:"
        )

        for error in errors:

            print(
                f"  - {error}"
            )


# ==========================================================
# ENTRY POINT
# ==========================================================

if __name__ == "__main__":

    main()