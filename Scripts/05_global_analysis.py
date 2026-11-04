#!/usr/bin/env python3

import os
import pandas as pd

from utils.statistics_utils import (
    protein_statistics,
    drug_statistics,
    mutation_statistics,
    affinity_summary,
    pocket_transition_matrix,
    distance_statistics,
    correlation_statistics,
)

from utils.plotting_utils import (
    plot_delta_distribution,
    plot_protein_barplot,
    plot_drug_barplot,
    plot_mutation_barplot,
    plot_transition_heatmap,
    plot_distance_distribution,
    plot_distance_vs_delta,
)

# ==========================================================
# INPUT
# ==========================================================

INPUT_CSV = "Analysis/03_compute_delta_scores/matched_pockets_scored.csv"

OUTPUT_DIR = "Results/Statistics"
FIGURE_DIR = "Results/Figures"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# ==========================================================
# LOAD
# ==========================================================

df = pd.read_csv(INPUT_CSV)

required = [
    "Protein",
    "Drug",
    "Mutation",
    "Distance",
    "WT_PocketRank",
    "Mut_PocketRank",
    "DeltaDockingScore",
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(f"Missing columns: {missing}")

# ==========================================================
# STATISTICS
# ==========================================================

print("Computing statistics...")

protein_stats = protein_statistics(df)
drug_stats = drug_statistics(df)
mutation_stats = mutation_statistics(df)

affinity_stats = affinity_summary(df)
transition_matrix = pocket_transition_matrix(df)

distance_stats = distance_statistics(df)
correlation_stats = correlation_statistics(df)

# ==========================================================
# SAVE TABLES
# ==========================================================

protein_stats.to_csv(
    os.path.join(OUTPUT_DIR, "protein_statistics.csv"),
    index=False,
)

drug_stats.to_csv(
    os.path.join(OUTPUT_DIR, "drug_statistics.csv"),
    index=False,
)

mutation_stats.to_csv(
    os.path.join(OUTPUT_DIR, "mutation_statistics.csv"),
    index=False,
)

affinity_stats.to_csv(
    os.path.join(OUTPUT_DIR, "affinity_summary.csv"),
    index=False,
)

transition_matrix.to_csv(
    os.path.join(OUTPUT_DIR, "pocket_transition_matrix.csv"),
)

distance_stats.to_csv(
    os.path.join(OUTPUT_DIR, "distance_statistics.csv"),
)

correlation_stats.to_csv(
    os.path.join(OUTPUT_DIR, "correlation_statistics.csv"),
    index=False,
)

# ==========================================================
# EXCEL WORKBOOK
# ==========================================================

with pd.ExcelWriter(
    os.path.join(OUTPUT_DIR, "analysis_summary.xlsx")
) as writer:

    protein_stats.to_excel(
        writer,
        sheet_name="Proteins",
        index=False,
    )

    drug_stats.to_excel(
        writer,
        sheet_name="Drugs",
        index=False,
    )

    mutation_stats.to_excel(
        writer,
        sheet_name="Mutations",
        index=False,
    )

    affinity_stats.to_excel(
        writer,
        sheet_name="Affinity",
        index=False,
    )

    transition_matrix.to_excel(
        writer,
        sheet_name="PocketTransitions",
    )

    distance_stats.to_excel(
        writer,
        sheet_name="Distance",
    )

    correlation_stats.to_excel(
        writer,
        sheet_name="Correlation",
        index=False,
    )

# ==========================================================
# FIGURES
# ==========================================================

print("Generating figures...")

plot_delta_distribution(df, FIGURE_DIR)

plot_protein_barplot(
    protein_stats,
    FIGURE_DIR,
)

plot_drug_barplot(
    drug_stats,
    FIGURE_DIR,
)

plot_mutation_barplot(
    mutation_stats,
    FIGURE_DIR,
)

plot_transition_heatmap(
    transition_matrix,
    FIGURE_DIR,
)

plot_distance_distribution(
    df,
    FIGURE_DIR,
)

plot_distance_vs_delta(
    df,
    FIGURE_DIR,
)

print()

print("Done!")

print()

print(f"Statistics -> {OUTPUT_DIR}")
print(f"Figures    -> {FIGURE_DIR}")