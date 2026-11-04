#!/usr/bin/env python3
"""
02_match_pockets.py

Matches WT and mutant pockets using structural alignment.

Inputs
------
analysis/alignment_matrices.pkl
top_3_pockets.json
pose_data.json

Outputs
-------
analysis/matched_pockets.csv
analysis/distance_matrix.csv
"""

from pathlib import Path
import pandas as pd
import numpy as np

from utils.pocket_matching import (
    alignments,
    pocket_json,
    pose_json,
    get_wildtype_key,
    split_ligand_name,
    best_docking_scores,
    build_pocket_list,
    transform_centers,
    build_distance_matrix,
    match_pockets,
)

###############################################################################
# OUTPUT
###############################################################################

OUTPUT_DIR = Path("Analysis/02_matched_pockets")
OUTPUT_DIR.mkdir(exist_ok=True)

MATCHED_OUTPUT = OUTPUT_DIR / "matched_pockets.csv"
DISTANCE_OUTPUT = OUTPUT_DIR / "distance_matrix.csv"

###############################################################################

matched_rows = []
distance_rows = []

###############################################################################
# Iterate over every Protein-Drug pair
###############################################################################

for ligand_name in sorted(pocket_json):

    protein, drug = split_ligand_name(ligand_name)

    print(f"\nProcessing {protein} - {drug}")

    pocket_entry = pocket_json[ligand_name]
    pose_entry = pose_json[ligand_name]

    ###########################################################################
    # Wild-type
    ###########################################################################

    wt_key = get_wildtype_key(pocket_entry)

    wt_structure = next(iter(pocket_entry[wt_key]))

    wt_docking = best_docking_scores(

        pose_entry[wt_key][wt_structure]

    )

    print("\n========================")
    print(ligand_name)
    print(wt_key)
    print(wt_structure)

    wt_pockets = build_pocket_list(

        pocket_entry[wt_key][wt_structure],

        wt_docking,

    )

    ###########################################################################
    # Mutants
    ###########################################################################

    for mutation in pocket_entry:

        if mutation == wt_key:
            continue

        alignment_key = (protein, mutation)

        if alignment_key not in alignments:

            print(
                f"Alignment missing for {protein} {mutation}"
            )

            continue

        rotation = alignments[alignment_key]["rotation"]
        translation = alignments[alignment_key]["translation"]

        mut_structure = next(

            iter(pocket_entry[mutation])

        )

        mut_docking = best_docking_scores(

            pose_entry[mutation][mut_structure]

        )

        mut_pockets = build_pocket_list(

            pocket_entry[mutation][mut_structure],

            mut_docking,

        )

        #######################################################################
        # Transform mutant centers
        #######################################################################

        mut_pockets = transform_centers(

            mut_pockets,

            rotation,

            translation,

        )

        #######################################################################
        # Distance matrix
        #######################################################################

        D = build_distance_matrix(

            wt_pockets,

            mut_pockets,

        )

        #######################################################################
        # Save complete matrix
        #######################################################################

        for i, wt in enumerate(wt_pockets):

            for j, mut in enumerate(mut_pockets):

                distance_rows.append(

                    {

                        "Protein": protein,

                        "Drug": drug,

                        "Mutation": mutation,

                        "WT_PocketRank": wt.rank,

                        "WT_PocketID": wt.pocket_id,

                        "Mut_PocketRank": mut.rank,

                        "Mut_PocketID": mut.pocket_id,

                        "Distance": float(D[i, j]),

                    }

                )

        #######################################################################
        # Hungarian assignment
        #######################################################################

        matches = match_pockets(D)

        #######################################################################
        # Final dataframe
        #######################################################################

        for wt_idx, mut_idx, dist in matches:

            wt = wt_pockets[wt_idx]
            mut = mut_pockets[mut_idx]

            matched_rows.append(

                {

                    "Protein": protein,

                    "Drug": drug,

                    "Mutation": mutation,

                    "WT_PocketRank": wt.rank,

                    "WT_PocketID": wt.pocket_id,

                    "Mut_PocketRank": mut.rank,

                    "Mut_PocketID": mut.pocket_id,

                    "Distance": dist,

                    "WT_PocketScore": wt.pocket_score,

                    "Mut_PocketScore": mut.pocket_score,

                    "WT_DockingScore": wt.docking_score,

                    "Mut_DockingScore": mut.docking_score,

                    "WT_X": wt.center[0],
                    "WT_Y": wt.center[1],
                    "WT_Z": wt.center[2],

                    "Mut_X": mut.center[0],
                    "Mut_Y": mut.center[1],
                    "Mut_Z": mut.center[2],

                }

            )

###############################################################################
# Save
###############################################################################

matched_df = pd.DataFrame(matched_rows)


dist = matched_df["Distance"]

print("\nPocket Matching Statistics")
print("---------------------------")
print(f"Number of matches: {len(dist)}")
print(f"Mean distance:     {dist.mean():.2f} Å")
print(f"Median distance:   {dist.median():.2f} Å")
print(f"Std distance:      {dist.std():.2f} Å")
print(f"Min distance:      {dist.min():.2f} Å")
print(f"Max distance:      {dist.max():.2f} Å")

print(f"≤1 Å : {(dist <= 1).mean()*100:.2f}%")
print(f"≤2 Å : {(dist <= 2).mean()*100:.2f}%")
print(f"≤3 Å : {(dist <= 3).mean()*100:.2f}%")
print(f"≤5 Å : {(dist <= 5).mean()*100:.2f}%")
print(f"≤5 Å : {(dist <= 6).mean()*100:.2f}%")

distance_df = pd.DataFrame(distance_rows)

matched_df.sort_values(

    [

        "Protein",

        "Drug",

        "Mutation",

        "WT_PocketRank",

    ],

    inplace=True,

)

distance_df.sort_values(

    [

        "Protein",

        "Drug",

        "Mutation",

        "WT_PocketRank",

        "Mut_PocketRank",

    ],

    inplace=True,

)

matched_df.to_csv(

    MATCHED_OUTPUT,

    index=False,

)

distance_df.to_csv(

    DISTANCE_OUTPUT,

    index=False,

)

print("\n========================================")
print("Finished.")
print(f"Matched pockets : {len(matched_df)}")
print(f"Distance matrix : {len(distance_df)}")
print("========================================")