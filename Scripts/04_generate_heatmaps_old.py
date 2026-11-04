#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.colors import TwoSlopeNorm


# ==========================================================
# INPUT / OUTPUT
# ==========================================================

INPUT_CSV = "Analysis/03_compute_delta_scores/docking_summary.csv"

OUTPUT_DIR = "Results/Heatmaps"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==========================================================
# LOAD DATA
# ==========================================================

df = pd.read_csv(INPUT_CSV)


# ==========================================================
# GLOBAL COLOR SCALE
# ==========================================================

delta_min = df["DeltaDockingScore"].min()
delta_max = df["DeltaDockingScore"].max()

print(f"Global range: {delta_min:.2f} to {delta_max:.2f}")

norm = TwoSlopeNorm(
    vmin=delta_min,
    vcenter=0,
    vmax=delta_max,
)


# ==========================================================
# ONE HEATMAP PER PROTEIN
# ==========================================================

for protein in sorted(df["Protein"].unique()):

    protein_df = df[df["Protein"] == protein]

    matrix = protein_df.pivot(
        index="Mutation",
        columns="Drug",
        values="DeltaDockingScore",
    )

    matrix = matrix.sort_index()

    matrix = matrix[
        sorted(matrix.columns)
    ]

    fig_width = max(8, len(matrix.columns) * 0.5)
    fig_height = max(5, len(matrix.index) * 0.45)

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height)
    )

    im = ax.imshow(
        matrix.values,
        cmap="coolwarm",
        norm=norm,
        aspect="auto",
    )

    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(
        matrix.columns,
        rotation=90,
        fontsize=8,
    )

    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(
        matrix.index,
        fontsize=8,
    )

    ax.set_title(
        protein,
        fontsize=16,
        weight="bold",
    )

    # ------------------------------------------------------
    # Annotate every cell
    # ------------------------------------------------------

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):

            value = matrix.iloc[i, j]

            if pd.notna(value):

                ax.text(
                    j,
                    i,
                    f"{value:.1f}",
                    ha="center",
                    va="center",
                    fontsize=6,
                )

    cbar = plt.colorbar(im)

    cbar.set_label(
        "Δ Docking Score (Mutant - WT)",
        fontsize=10,
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"{protein}_heatmap.png",
        ),
        dpi=300,
        bbox_inches="tight",
    )

    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"{protein}_heatmap.pdf",
        ),
        bbox_inches="tight",
    )

    plt.close()

print("\nDone!")