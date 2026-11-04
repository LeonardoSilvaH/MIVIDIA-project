#!/usr/bin/env python3

"""
01_align_structures.py

Align all mutant AlphaFold models against their corresponding WT model.

Directory expected:

models/
├── MYC/
│   ├── wild-type/
│   │   └── MYC.pdb
│   ├── pP74A/
│   │   └── ranked_0.pdb
│   ├── pT58I/
│   │   └── ranked_0.pdb
│   └── ...
├── KRAS/
│   ├── wild-type/
│   │   └── KRAS.pdb
│   └── pG12D/
│       └── ranked_0.pdb

Outputs

analysis/
    alignment_results.csv
    alignment_matrices.pkl
"""

from pathlib import Path
import pickle
import pandas as pd
from tqdm import tqdm

from utils.alignment import (
    parse_mutation,
    superimpose_structures,
)

###############################################################################
# USER SETTINGS
###############################################################################

MODEL_DIR = Path("Data/models")

OUTPUT_DIR = Path("Analysis/01_align_structures")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_OUTPUT = OUTPUT_DIR / "alignment_results.csv"

PKL_OUTPUT = OUTPUT_DIR / "alignment_matrices.pkl"

###############################################################################
# MAIN
###############################################################################

results = []

matrices = {}

protein_folders = sorted(

    d

    for d in MODEL_DIR.iterdir()

    if d.is_dir()

)

print(f"\nFound {len(protein_folders)} proteins.\n")

for protein_dir in tqdm(protein_folders):

    protein = protein_dir.name

    wt_folder = protein_dir / "wild-type"

    if not wt_folder.exists():

        print(f"Skipping {protein}: WT folder not found.")

        continue

    wt_files = list(wt_folder.glob("*.pdb"))

    if len(wt_files) == 0:

        print(f"Skipping {protein}: WT pdb not found.")

        continue

    wt_pdb = wt_files[0]

    mutation_folders = [

        d

        for d in protein_dir.iterdir()

        if d.is_dir()

        and d.name != "wild-type"

    ]

    for mutation_folder in mutation_folders:

        mutation = mutation_folder.name

        # Remove AlphaFold prefix
        if mutation.startswith("p"):
            mutation = mutation[1:]

        # Ignore wild type (extra safety)
        if mutation.lower() == "wild-type":
            continue

        # If the folder doesn't already contain the protein name,
        # prepend it to match the CoBDock convention.
        if not mutation.startswith(f"{protein}-"):
            mutation = f"{protein}-{mutation}"

        if mutation is None:

            continue

        pdb_file = mutation_folder / "ranked_0.pdb"

        if not pdb_file.exists():

            print(f"{protein} {mutation}: ranked_0.pdb not found.")

            continue

        try:

            alignment = superimpose_structures(

                wt_pdb=str(wt_pdb),

                mut_pdb=str(pdb_file),

                mutation=mutation,

            )

        except Exception as e:

            print(

                f"Failed: {protein} {mutation}"

            )

            print(e)

            continue

        results.append(

            {

                "Protein": protein,

                "Mutation": mutation,

                "WT_PDB": str(wt_pdb),

                "Mutant_PDB": str(pdb_file),

                "RMSD": alignment["rmsd"],

                "AlignedResidues": alignment["aligned_residues"],

            }

        )

        matrices[(protein, mutation)] = {

            "rotation": alignment["rotation"],

            "translation": alignment["translation"],

            "rmsd": alignment["rmsd"],

            "aligned_residues": alignment["aligned_residues"],

            "residue_numbers": alignment["residue_numbers"],

        }

###############################################################################
# SAVE RESULTS
###############################################################################

print("\n===== Alignment keys =====")
for key in list(alignment.keys())[:20]:
    print(key)

df = pd.DataFrame(results)

df.sort_values(

    ["Protein", "Mutation"],

    inplace=True,

)

df.to_csv(

    CSV_OUTPUT,

    index=False,

)

with open(PKL_OUTPUT, "wb") as f:

    pickle.dump(

        matrices,

        f,

        protocol=pickle.HIGHEST_PROTOCOL,

    )

print()

print("=" * 60)

print("Alignment finished.")

print(f"Proteins processed : {df['Protein'].nunique()}")

print(f"Mutations processed: {len(df)}")

print()

print(f"CSV saved to : {CSV_OUTPUT}")

print(f"Matrices saved to : {PKL_OUTPUT}")

print("=" * 60)