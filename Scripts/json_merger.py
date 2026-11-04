#!/usr/bin/env python3

import shutil
import json
from pathlib import Path

# ==========================================================
# INPUTS
# ==========================================================

RUN1 = Path("Data/results_KRAS/pocket_locations/")
RUN2 = Path("Data/results_semifinal/pocket_locations/")

OUTPUT = Path("Data/results_merged/")
OUTPUT.mkdir(exist_ok=True)

# ==========================================================
# FUNCTIONS
# ==========================================================

def merge_json(file1, file2, output_file):

    with open(file1) as f:
        data1 = json.load(f)

    with open(file2) as f:
        data2 = json.load(f)

    merged = {}
    duplicates = []

    # Add first run
    for key, value in data1.items():
        merged[key] = value

    # Add second run
    for key, value in data2.items():

        if key in merged:
            duplicates.append(key)
            print(f"Duplicate found: {key}")
            continue

        merged[key] = value

    with open(output_file, "w") as f:
        json.dump(merged, f, indent=4)

    return merged, duplicates

def merge_local_docking(run1, run2, output):

    src1 = run1 / "local_docking"
    src2 = run2 / "local_docking"
    dst = output / "local_docking"

    dst.mkdir(exist_ok=True)

    copied = 0
    duplicates = []

    for source in (src1, src2):

        if not source.exists():
            continue

        for pair_dir in source.iterdir():

            if not pair_dir.is_dir():
                continue

            destination = dst / pair_dir.name

            if destination.exists():
                duplicates.append(pair_dir.name)
                print(f"Duplicate local_docking folder: {pair_dir.name}")
                continue

            shutil.copytree(pair_dir, destination)
            copied += 1

    return copied, duplicates

# ==========================================================
# MERGE
# ==========================================================

print("Merging top_3_pockets.json...")

top3, dup_top = merge_json(
    RUN1 / "top_3_pockets.json",
    RUN2 / "top_3_pockets.json",
    OUTPUT / "top_3_pockets.json",
)

print()

print("Merging pose_data.json...")

poses, dup_pose = merge_json(
    RUN1 / "pose_data.json",
    RUN2 / "pose_data.json",
    OUTPUT / "pose_data.json",
)

print()

print("Merging local_docking...")

copied, dup_dirs = merge_local_docking(
    RUN1,
    RUN2,
    OUTPUT,
)

# ==========================================================
# VALIDATION
# ==========================================================

top_keys = set(top3.keys())
pose_keys = set(poses.keys())

missing_in_pose = top_keys - pose_keys
missing_in_top = pose_keys - top_keys

print("\n========================")

print(f"Total Protein-Drug pairs : {len(top_keys)}")

print(f"Duplicates skipped (top_3): {len(dup_top)}")
print(f"Duplicates skipped (poses): {len(dup_pose)}")

print()

print(f"Protein-drug folders copied : {copied}")
print(f"Duplicate folders skipped   : {len(dup_dirs)}")

if missing_in_pose:
    print("Present in top_3_pockets but missing in pose_data:")
    for k in sorted(missing_in_pose):
        print("  ", k)

if missing_in_top:
    print("Present in pose_data but missing in top_3_pockets:")
    for k in sorted(missing_in_top):
        print("  ", k)

if not missing_in_pose and not missing_in_top:
    print("✓ Validation successful.")
    print("✓ Both JSON files contain identical Protein-Drug pairs.")

print(f"\nMerged files written to:\n{OUTPUT}")