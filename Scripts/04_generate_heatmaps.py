#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.colors import TwoSlopeNorm


# ==========================================================
# INPUT / OUTPUT
# ==========================================================

INPUT_CSV = (
    "Analysis/03_compute_delta_scores/docking_summary.csv"
)

OUTPUT_DIR = (
    "Results/Heatmaps"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


# ==========================================================
# LOAD DATA
# ==========================================================

print("Loading data...")

df = pd.read_csv(
    INPUT_CSV
)

print(
    f"Loaded {len(df)} rows."
)


# ==========================================================
# SCORE
# ==========================================================

SCORE_COLUMN = "DeltaDockingScore"

SCORE_LABEL = (
    "Δ Docking Score (Mutant − WT)"
)

SCORE_TYPE = "PLANTS"

if SCORE_COLUMN not in df.columns:

    raise KeyError(
        f"\nColumn '{SCORE_COLUMN}' not found in {INPUT_CSV}.\n\n"
        f"Available columns:\n{list(df.columns)}"
    )


# ==========================================================
# GLOBAL COLOR SCALE
# ==========================================================

delta_min = df[
    SCORE_COLUMN
].min()

delta_max = df[
    SCORE_COLUMN
].max()

print(
    f"Global range: "
    f"{delta_min:.2f} to "
    f"{delta_max:.2f}"
)


# ----------------------------------------------------------
# Handle case where all values are on one side of zero
# ----------------------------------------------------------

if delta_min < 0 and delta_max > 0:

    norm = TwoSlopeNorm(
        vmin=delta_min,
        vcenter=0,
        vmax=delta_max,
    )

else:

    print(
        "WARNING: Score range does not cross zero."
    )

    # Expand the range slightly so that the heatmap
    # can still be generated.

    maximum = max(
        abs(delta_min),
        abs(delta_max),
        1
    )

    norm = TwoSlopeNorm(
        vmin=-maximum,
        vcenter=0,
        vmax=maximum,
    )


# ==========================================================
# ONE HEATMAP PER PROTEIN
# ==========================================================

for protein in sorted(
    df["Protein"].dropna().unique()
):

    print(
        f"\nGenerating heatmap for {protein}..."
    )

    protein_df = df[
        df["Protein"] == protein
    ].copy()


    # ======================================================
    # CREATE MUTATION × DRUG MATRIX
    # ======================================================

    matrix = protein_df.pivot_table(
        index="Mutation",
        columns="Drug",
        values=SCORE_COLUMN,
        aggfunc="mean",
    )


    # ------------------------------------------------------
    # Sort rows and columns
    # ------------------------------------------------------

    matrix = matrix.sort_index()

    matrix = matrix[
        sorted(
            matrix.columns
        )
    ]


    # ======================================================
    # FIGURE SIZE
    # ======================================================

    fig_width = max(
        8,
        len(matrix.columns) * 0.5
    )

    fig_height = max(
        5,
        len(matrix.index) * 0.45
    )


    fig, ax = plt.subplots(
        figsize=(
            fig_width,
            fig_height
        )
    )


    # ======================================================
    # HEATMAP
    # ======================================================

    im = ax.imshow(
        matrix.values,
        cmap="coolwarm",
        norm=norm,
        aspect="auto",
    )


    # ======================================================
    # X AXIS
    # ======================================================

    ax.set_xticks(
        range(
            len(
                matrix.columns
            )
        )
    )

    ax.set_xticklabels(
        matrix.columns,
        rotation=90,
        fontsize=8,
    )


    # ======================================================
    # Y AXIS
    # ======================================================

    ax.set_yticks(
        range(
            len(
                matrix.index
            )
        )
    )

    ax.set_yticklabels(
        matrix.index,
        fontsize=8,
    )


    # ======================================================
    # TITLE
    # ======================================================

    ax.set_title(
        f"{protein}",
        fontsize=16,
        weight="bold",
    )


    # ======================================================
    # ANNOTATE CELLS
    # ======================================================

    for i in range(
        matrix.shape[0]
    ):

        for j in range(
            matrix.shape[1]
        ):

            value = matrix.iloc[
                i,
                j
            ]

            if pd.notna(
                value
            ):

                ax.text(
                    j,
                    i,
                    f"{value:.1f}",
                    ha="center",
                    va="center",
                    fontsize=6,
                )


    # ======================================================
    # COLORBAR
    # ======================================================

    cbar = plt.colorbar(
        im,
        ax=ax
    )

    cbar.set_label(
        SCORE_LABEL,
        fontsize=10,
    )


    # ======================================================
    # LAYOUT
    # ======================================================

    plt.tight_layout()


    # ======================================================
    # SAVE PNG
    # ======================================================

    png_path = os.path.join(
        OUTPUT_DIR,
        f"{protein}_{SCORE_TYPE}_heatmap.png"
    )

    plt.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )


    # ======================================================
    # SAVE PDF
    # ======================================================

    pdf_path = os.path.join(
        OUTPUT_DIR,
        f"{protein}_{SCORE_TYPE}_heatmap.pdf"
    )

    plt.savefig(
        pdf_path,
        bbox_inches="tight",
    )


    plt.close()


    print(
        f"  PNG: {png_path}"
    )

    print(
        f"  PDF: {pdf_path}"
    )


# ==========================================================
# DONE
# ==========================================================

print(
    "\nDone!"
)

print(
    f"Generated {len(df['Protein'].unique())} "
    f"protein heatmap(s)."
)