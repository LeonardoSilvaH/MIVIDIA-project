import pandas as pd
from scipy.stats import pearsonr, spearmanr


# ==========================================================
# Internal score helper
# ==========================================================

def prepare_delta(df):
    """
    Return a copy of the dataframe with a standardized
    internal column named DeltaScore.
    """

    df = df.copy()

    df["DeltaScore"] = df[
        "DeltaDockingScore"
    ]

    return df


# ==========================================================
# Protein statistics
# ==========================================================

def protein_statistics(df):

    df = prepare_delta(df)

    return (
        df.groupby("Protein")
        .agg(
            NDrugs=("Drug", "nunique"),
            NMutations=("Mutation", "nunique"),

            MeanDelta=(
                "DeltaScore",
                "mean"
            ),

            MeanAbsDelta=(
                "DeltaScore",
                lambda x: x.abs().mean()
            ),

            Median=(
                "DeltaScore",
                "median"
            ),

            Std=(
                "DeltaScore",
                "std"
            ),

            MaxAbsDelta=(
                "DeltaScore",
                lambda x: x.abs().max()
            ),
        )
        .reset_index()
        .sort_values(
            "MeanAbsDelta",
            ascending=False
        )
    )


# ==========================================================
# Drug statistics
# ==========================================================

def drug_statistics(df):

    df = prepare_delta(df)

    return (
        df.groupby("Drug")
        .agg(
            NProteins=(
                "Protein",
                "nunique"
            ),

            NMutations=(
                "Mutation",
                "nunique"
            ),

            MeanDelta=(
                "DeltaScore",
                "mean"
            ),

            MeanAbsDelta=(
                "DeltaScore",
                lambda x: x.abs().mean()
            ),

            Std=(
                "DeltaScore",
                "std"
            ),

            MaxAbsDelta=(
                "DeltaScore",
                lambda x: x.abs().max()
            ),
        )
        .reset_index()
        .sort_values(
            "MeanAbsDelta",
            ascending=False
        )
    )


# ==========================================================
# Mutation statistics
# ==========================================================

def mutation_statistics(df):

    df = prepare_delta(df)

    return (
        df.groupby(
            [
                "Protein",
                "Mutation"
            ]
        )
        .agg(
            NDrugs=(
                "Drug",
                "nunique"
            ),

            MeanDelta=(
                "DeltaScore",
                "mean"
            ),

            MeanAbsDelta=(
                "DeltaScore",
                lambda x: x.abs().mean()
            ),

            Std=(
                "DeltaScore",
                "std"
            ),

            MaxAbsDelta=(
                "DeltaScore",
                lambda x: x.abs().max()
            ),
        )
        .reset_index()
        .sort_values(
            "MeanAbsDelta",
            ascending=False
        )
    )


# ==========================================================
# Affinity summary
# ==========================================================

def classify_affinity(delta):

    if delta <= -1:

        return "Stronger"

    elif delta >= 1:

        return "Weaker"

    return "Neutral"


def affinity_summary(df):

    df = prepare_delta(df)

    tmp = df.copy()

    tmp["AffinityChange"] = (
        tmp["DeltaScore"]
        .apply(classify_affinity)
    )

    summary = (
        tmp["AffinityChange"]
        .value_counts()
        .rename_axis("Class")
        .reset_index(
            name="Count"
        )
    )

    summary["Percentage"] = (
        summary["Count"]
        /
        summary["Count"].sum()
        *
        100
    )

    return summary


# ==========================================================
# Pocket transitions
# ==========================================================

def pocket_transition_matrix(df):

    return pd.crosstab(
        df["WT_PocketRank"],
        df["Mut_PocketRank"],
    )


# ==========================================================
# Distance statistics
# ==========================================================

def distance_statistics(df):

    stats = (
        df["Distance"]
        .describe()
        .to_frame()
        .T
    )

    return stats


# ==========================================================
# Correlations
# ==========================================================

def correlation_statistics(df):

    df = prepare_delta(df)

    # ------------------------------------------------------
    # Remove missing values
    # ------------------------------------------------------

    correlation_df = df[
        [
            "Distance",
            "DeltaScore"
        ]
    ].dropna()

    if len(correlation_df) < 2:

        return pd.DataFrame(
            {
                "Metric": [
                    "Pearson",
                    "Spearman"
                ],

                "Correlation": [
                    float("nan"),
                    float("nan")
                ],

                "PValue": [
                    float("nan"),
                    float("nan")
                ]
            }
        )

    abs_delta = (
        correlation_df[
            "DeltaScore"
        ].abs()
    )

    # ------------------------------------------------------
    # Pearson
    # ------------------------------------------------------

    pearson_r, pearson_p = pearsonr(
        correlation_df["Distance"],
        abs_delta,
    )

    # ------------------------------------------------------
    # Spearman
    # ------------------------------------------------------

    spearman_rho, spearman_p = spearmanr(
        correlation_df["Distance"],
        abs_delta,
    )

    return pd.DataFrame(
        {
            "Metric": [
                "Pearson",
                "Spearman",
            ],

            "Correlation": [
                pearson_r,
                spearman_rho,
            ],

            "PValue": [
                pearson_p,
                spearman_p,
            ],
        }
    )