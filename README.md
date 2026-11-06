# MIVIDIA project

Mutational Impact of Variations In Drug Interaction Analysis.

Compares how point mutations change drug binding by docking each drug against
wild-type and mutant AlphaFold models (CoBDock / PLANTS) and comparing the
docking scores **within the same structural pocket**.

Continued from the
[KRAS article repository](https://github.com/LeonardoSilvaH/MIVIDIA-Mutational-Impact-of-Variations-In-Drug-Interaction-Analysis),
which is kept unchanged as published.

## Pipeline

Run every script from the repository root.

| Step | Script | What it does |
|---|---|---|
| 1 | `Scripts/01_align_structures.py` | Superimposes each mutant model onto its wild type (Cα, chain A) |
| 2 | `Scripts/02_match_pockets.py` | Matches WT and mutant pockets by residue overlap; labels each WT pocket Matched / NotDocked / Lost / Unresolved |
| 3 | `Scripts/03_compute_delta_scores.py` | Δ PLANTS score (mutant − WT) for Matched pockets only |
| 4 | `Scripts/04_generate_heatmaps.py` | Mutation × drug heatmap per protein |
| 5 | `Scripts/05_global_analysis.py` | Summary statistics and figures |

## Paths

Every script resolves its paths from the repository itself, so they can be
launched from any working directory:

```bash
python /hd/leo/MIVIDIA-project/Scripts/02_match_pockets.py
```

`Analysis/` and `Results/` are always written inside the repository.

## Inputs

- `Data/models/<PROTEIN>/wild-type/*.pdb` and
  `Data/models/<PROTEIN>/p<MUT>/ranked_0.pdb` (in this repository)
- CoBDock output, read from `$MIVIDIA_RESULTS` (default
  `/hd/leo/MIVIDIA/Data/resultados_mestrado`):
  - `pocket_locations/top_3_pockets.json`, `pocket_locations/pose_data.json`
  - `docking/fpocket_data.json`, `docking/p2rank_data.json` (pocket residues)

To analyse another results folder, point the variable at it instead of
editing the scripts:

```bash
export MIVIDIA_RESULTS=/path/to/cobdock/results
```

## Setup

```bash
conda create -p <env-path> -c conda-forge --override-channels python=3.11 pandas numpy scipy matplotlib biopython tqdm openpyxl
```
