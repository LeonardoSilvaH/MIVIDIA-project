#!/usr/bin/env python3

import json
import pandas as pd

# ==========================================================
# INPUT / OUTPUT
# ==========================================================

POSE_JSON = "Data/results_merged/pocket_locations/pose_data.json"
OUTPUT_CSV = "Analysis/     q1.5_Docking_csv/all_docking_poses.csv"

# ==========================================================
# LOAD JSON
# ==========================================================

with open(POSE_JSON) as f:
    data = json.load(f)

rows = []

# ==========================================================
# PARSE
# ==========================================================

for protein_drug, mutations in data.items():

    # Split "BRAF_Gilteritinib"
    protein, drug = protein_drug.split("_", 1)

    for mutation_name, chains in mutations.items():

        mutation = mutation_name

        if "-" in mutation:
            mutation = mutation.split("-")[-1]

        for chain_name, pockets in chains.items():

            for pocket_id, poses in pockets.items():

                for pose_id, pose in poses.items():

                    row = {

                        "Protein": protein,
                        "Drug": drug,
                        "Mutation": mutation,
                        "MutationFolder": mutation_name,
                        "Chain": chain_name,

                        "PocketID": pocket_id,
                        "PoseID": pose_id,

                        "CenterX": pose.get("center_x"),
                        "CenterY": pose.get("center_y"),
                        "CenterZ": pose.get("center_z"),

                        "SizeX": pose.get("size_x"),
                        "SizeY": pose.get("size_y"),
                        "SizeZ": pose.get("size_z"),

                        "DockingScore": pose.get("total_score"),

                        "Score_RB_Pen": pose.get("score_rb_pen"),
                        "Score_Norm_HeavyAtoms": pose.get("score_norm_hevatoms"),
                        "Score_Norm_CorrectedHeavyAtoms": pose.get("score_norm_crt_hevatoms"),
                        "Score_Norm_Weight": pose.get("score_norm_weight"),
                        "Score_Norm_CorrectedWeight": pose.get("score_norm_crt_weight"),
                        "Score_RB_Pen_CorrectedHeavyAtoms": pose.get("score_rb_pen_norm_crt_hevatoms"),
                        "Score_Norm_Contact": pose.get("score_norm_contact"),

                    }

                    rows.append(row)

# ==========================================================
# SAVE
# ==========================================================

df = pd.DataFrame(rows)

df = df.sort_values(
    [
        "Protein",
        "Drug",
        "Mutation",
        "PocketID",
        "PoseID",
    ]
)

df.to_csv(OUTPUT_CSV, index=False)

print()
print(f"Extracted {len(df):,} docking poses.")
print(f"Saved to {OUTPUT_CSV}")