#!/usr/bin/env python3

import json
import pandas as pd
from tqdm import tqdm


# ==========================================================
# INPUT / OUTPUT
# ==========================================================

MATCHED_POCKETS_FILE = "Analysis/02_matched_pockets/matched_pockets.csv"
POSE_DATA_FILE = "Data/results_KRAS/pocket_locations/pose_data.json"

OUTPUT_FULL = "Analysis/03_compute_delta_scores/matched_pockets_scored.csv"
OUTPUT_SUMMARY = "Analysis/03_compute_delta_scores/docking_summary.csv"


# ==========================================================
# HELPERS
# ==========================================================

def load_json(filename):
    with open(filename, "r") as f:
        return json.load(f)


def get_chain_name(mutation_dict):
    """Returns the only chain available (usually *_A)."""
    return next(iter(mutation_dict.keys()))


def get_best_pose(pocket_dict):
    """
    Returns
    -------
    best_pose : str
    best_score : float
    """

    best_pose = None
    best_score = float("inf")

    for pose_id, pose_data in pocket_dict.items():

        score = pose_data["total_score"]

        if score < best_score:
            best_score = score
            best_pose = pose_id

    return best_pose, best_score


def get_pocket(
    pose_data,
    protein,
    drug,
    mutation,
    pocket_rank,
):
    """
    Returns the dictionary corresponding to one pocket.
    """

    complex_name = f"{protein}_{drug}"

    mutation_dict = pose_data[complex_name][mutation]

    chain = get_chain_name(mutation_dict)

    return mutation_dict[chain][str(int(pocket_rank))]


# ==========================================================
# PROCESS ONE MATCH
# ==========================================================

def process_match(row, pose_data):

    protein = row["Protein"]
    drug = row["Drug"]
    mutation = row["Mutation"]

    wt_mutation = f"{protein}-WILD-TYE"

    wt_rank = row["WT_PocketRank"]
    mut_rank = row["Mut_PocketRank"]

    wt_pocket = get_pocket(
        pose_data,
        protein,
        drug,
        wt_mutation,
        wt_rank,
    )

    mut_pocket = get_pocket(
        pose_data,
        protein,
        drug,
        mutation,
        mut_rank,
    )

    wt_pose, wt_score = get_best_pose(wt_pocket)
    mut_pose, mut_score = get_best_pose(mut_pocket)

    row["WT_BestPose"] = int(wt_pose)
    row["WT_BestDockingScore"] = wt_score

    row["Mut_BestPose"] = int(mut_pose)
    row["Mut_BestDockingScore"] = mut_score

    row["DeltaDockingScore"] = mut_score - wt_score

    row["AbsDeltaDockingScore"] = abs(row["DeltaDockingScore"])

    if row["DeltaDockingScore"] < 0:
        row["AffinityChange"] = "Stronger"

    elif row["DeltaDockingScore"] > 0:
        row["AffinityChange"] = "Weaker"

    else:
        row["AffinityChange"] = "Unchanged"

    return row


# ==========================================================
# MAIN
# ==========================================================

def main():

    print("Loading data...")

    matched_df = pd.read_csv(MATCHED_POCKETS_FILE)

    pose_data = load_json(POSE_DATA_FILE)

    processed_rows = []

    print(f"Processing {len(matched_df)} matched pockets...\n")

    for _, row in tqdm(
        matched_df.iterrows(),
        total=len(matched_df),
    ):

        try:

            processed_rows.append(
                process_match(
                    row.copy(),
                    pose_data,
                )
            )

        except Exception as e:

            print(
                f"Skipping "
                f"{row['Protein']} "
                f"{row['Mutation']} "
                f"{row['Drug']}"
            )

            print(e)

    scored_df = pd.DataFrame(processed_rows)

    scored_df = scored_df.sort_values(
        [
            "Protein",
            "Drug",
            "Mutation",
            "WT_PocketRank",
        ]
    )

    print("\nSaving complete table...")

    scored_df.to_csv(
        OUTPUT_FULL,
        index=False,
    )

    print("Generating summary table...")

    summary_df = (
        scored_df
        .sort_values(
            "Mut_BestDockingScore"
        )
        .groupby(
            [
                "Protein",
                "Drug",
                "Mutation",
            ],
            as_index=False,
        )
        .first()
    )

    summary_df.to_csv(
        OUTPUT_SUMMARY,
        index=False,
    )

    print("\nDone.")
    print(f"Saved: {OUTPUT_FULL}")
    print(f"Saved: {OUTPUT_SUMMARY}")


if __name__ == "__main__":
    main()