#!/usr/bin/env python3

"""
02b_matching_report.py

Report on the pocket matching, so MIN_RESIDUE_OVERLAP is chosen from
the data instead of assumed.

Answers four questions:

  1. Is the residue-overlap distribution bimodal, with the cutoff in an
     empty gap, or is it continuous across the cutoff?
  2. How many Lost pockets are near misses - remodelled by the mutation
     rather than genuinely absent?
  3. Does comparing an fpocket pocket against a p2rank pocket deflate
     the Jaccard overlap, and does the overlap coefficient recover it?
  4. Which matched pairs sit far apart after superposition, i.e. where
     the alignment disagrees with the residue evidence?

Read-only. Run after 02_match_pockets.py.
"""

from pathlib import Path
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

MATCHED_FILE = (
    PROJECT_ROOT /
    "Analysis" /
    "02_matched_pockets" /
    "02_matched_pockets.csv"
)

MIN_RESIDUE_OVERLAP = 0.5

CUTOFFS = (0.3, 0.4, 0.5, 0.6, 0.7)


def banner(number, title):

    print()
    print("=" * 74)
    print(f"{number}. {title}")
    print("=" * 74)


def histogram(values, width=48):
    """
    Text histogram over 0-1 in 0.05 bins.
    """

    values = values[~np.isnan(values)]

    if not len(values):

        print("  (no values)")

        return

    counts, edges = np.histogram(values, bins=20, range=(0.0, 1.0))

    peak = counts.max() or 1

    for count, low in zip(counts, edges[:-1]):

        bar = "#" * int(width * count / peak)

        mark = "  <== cutoff" if abs(low - MIN_RESIDUE_OVERLAP) < 1e-9 else ""

        print(f"  {low:.2f}-{low + 0.05:.2f} {count:7d} |{bar}{mark}")


df = pd.read_csv(MATCHED_FILE)

print(f"Read {len(df)} rows from {MATCHED_FILE}")

print()
print(df["MatchStatus"].value_counts().to_string())

# ---------------------------------------------------------------- 1

banner(1, "RESIDUE OVERLAP DISTRIBUTION")

for status in ("Matched", "NotDocked", "Lost"):

    subset = df.loc[df["MatchStatus"] == status, "ResidueOverlap"]

    print(f"\n  {status}  (n={len(subset)}, "
          f"{subset.notna().sum()} with an overlap value)")

    histogram(subset.to_numpy(dtype=float))

print("\n  If the Lost bars run right up to the cutoff, the cutoff is")
print("  splitting one continuous population rather than separating two.")

# ---------------------------------------------------------------- 2

banner(2, "NEAR MISSES AMONG Lost")

lost = df[df["MatchStatus"] == "Lost"].copy()

near = lost[lost["ResidueOverlap"].notna()]

print(f"  Lost pockets                 : {len(lost)}")
print(f"  with a best-overlap recorded : {len(near)}")

if len(near):

    for low in (0.3, 0.4, 0.45):

        n = (near["ResidueOverlap"] >= low).sum()

        print(f"    best overlap >= {low:.2f}     : {n} "
              f"({100 * n / len(lost):.1f}% of Lost)")

    close = near[
        (near["ResidueOverlap"] >= 0.3)
        & (near["Distance"] <= 3.0)
    ]

    print(f"\n  overlap >= 0.30 AND within 3 A : {len(close)}")
    print("  These are almost certainly the same site, remodelled.")

    if len(close):

        print("\n  most affected proteins:")

        for protein, n in close["Protein"].value_counts().head(8).items():

            print(f"    {protein:10s} {n}")

else:

    print("\n  No overlap values on Lost rows - rerun 02_match_pockets.py")
    print("  with the current version, which records them.")

# ---------------------------------------------------------------- 3

banner(3, "CUTOFF SENSITIVITY")

scored = df[df["ResidueOverlap"].notna()]

print(f"  Pairs with an overlap value: {len(scored)}\n")
print(f"  {'cutoff':>8}  {'would match':>12}  {'change':>8}")

base = None

for cutoff in CUTOFFS:

    n = (scored["ResidueOverlap"] >= cutoff).sum()

    if base is None:

        base = n

    marker = "  <== current" if cutoff == MIN_RESIDUE_OVERLAP else ""

    print(f"  {cutoff:8.2f}  {n:12d}  {n - base:+8d}{marker}")

# ---------------------------------------------------------------- 4

banner(4, "DETECTOR MISMATCH")

if "SameProgram" in df.columns and "OverlapCoefficient" in df.columns:

    pairs = df[df["ResidueOverlap"].notna() & df["SameProgram"].notna()]

    for same in (True, False):

        subset = pairs[pairs["SameProgram"].astype(str) == str(same)]

        if not len(subset):

            continue

        label = "same detector" if same else "fpocket vs p2rank"

        print(f"\n  {label}  (n={len(subset)})")
        print(f"    Jaccard overlap     : mean "
              f"{subset['ResidueOverlap'].mean():.3f}  "
              f"median {subset['ResidueOverlap'].median():.3f}")
        print(f"    overlap coefficient : mean "
              f"{subset['OverlapCoefficient'].mean():.3f}  "
              f"median {subset['OverlapCoefficient'].median():.3f}")

        gain = (
            subset["OverlapCoefficient"] - subset["ResidueOverlap"]
        ).mean()

        print(f"    coefficient - Jaccard: {gain:+.3f}")

    print("\n  A large gap on the mixed-detector rows means the Jaccard")
    print("  overlap is being penalised for detector disagreement rather")
    print("  than for any difference between the structures.")

else:

    print("  Columns missing - rerun 02_match_pockets.py.")

# ---------------------------------------------------------------- 5

banner(5, "ALIGNMENT DISAGREEMENT (Matched pairs)")

matched = df[
    (df["MatchStatus"] == "Matched")
    & df["Distance"].notna()
]

if len(matched):

    far = matched[matched["Distance"] > 5.0]

    print(f"  Matched pairs            : {len(matched)}")
    print(f"  further than 5 A apart   : {len(far)} "
          f"({100 * len(far) / len(matched):.1f}%)")
    print(f"  median distance          : {matched['Distance'].median():.2f} A")
    print(f"  maximum distance         : {matched['Distance'].max():.2f} A")

    if len(far):

        print("\n  by protein (residues agree, superposition does not):")

        summary = (
            far.groupby("Protein")["Distance"]
            .agg(["size", "median", "max"])
            .sort_values("size", ascending=False)
            .head(10)
        )

        print(f"\n    {'protein':10s} {'n':>5s} {'median':>8s} {'max':>8s}")

        for protein, row in summary.iterrows():

            print(f"    {protein:10s} {int(row['size']):5d} "
                  f"{row['median']:8.2f} {row['max']:8.2f}")

        print("\n  These rely on residue identity, which is why the pipeline")
        print("  does not match on distance. Worth confirming the alignment")
        print("  RMSD for these proteins in 01_align_structures.")


# ---------------------------------------------------------------- 6

banner(6, "EVIDENCE STRENGTH (how many residues each match rests on)")

if "SharedResidues" in df.columns:

    ev = df[
        (df["MatchStatus"] == "Matched")
        & df["SharedResidues"].notna()
    ].copy()

    ev["SharedResidues"] = ev["SharedResidues"].astype(int)

    print(f"  Matched pairs with a shared-residue count: {len(ev)}\n")

    print(f"  {'shared':>8}  {'pairs':>8}  {'cumulative %':>13}")

    counts = ev["SharedResidues"].value_counts().sort_index()

    cumulative = 0

    for shared, n in counts.items():

        if shared > 12:

            break

        cumulative += n

        print(f"  {shared:8d}  {n:8d}  {100 * cumulative / len(ev):12.1f}%")

    big = (ev["SharedResidues"] > 12).sum()

    print(f"  {'>12':>8}  {big:8d}  {100.0:12.1f}%")

    print("\n  A pair sharing 2 residues out of 3 in each pocket reaches a")
    print("  Jaccard of 0.50 and is accepted on the same footing as a pair")
    print("  sharing 33 residues out of 50. These are not equal evidence.")

    print("\n  Effect of an absolute floor on shared residues:\n")

    print(f"  {'floor':>7}  {'kept':>8}  {'dropped':>8}  {'% dropped':>10}")

    for floor in (3, 4, 5, 6, 8):

        kept = (ev["SharedResidues"] >= floor).sum()

        dropped = len(ev) - kept

        print(f"  {floor:7d}  {kept:8d}  {dropped:8d}  "
              f"{100 * dropped / len(ev):9.1f}%")

    print("\n  Smallest pocket in each matched pair:\n")

    sizes = df.loc[ev.index, ["WT_NResidues", "Mut_NResidues"]].min(axis=1)

    for low, high in ((0, 4), (4, 6), (6, 9), (9, 15), (15, 10000)):

        n = ((sizes >= low) & (sizes < high)).sum()

        label = f"{low}-{high - 1}" if high < 10000 else f"{low}+"

        print(f"    {label:>8} residues : {n:6d} pairs "
              f"({100 * n / len(ev):5.1f}%)")

else:

    print("  SharedResidues column missing - rerun 02_match_pockets.py.")
