#!/usr/bin/env python3
"""
02_match_pockets.py

Match WT and mutant pockets after structural alignment.

Author:
    Leo Silva + ChatGPT

Inputs
------
analysis/alignment_matrices.pkl
pose_data.json
top_3_pockets.json

Outputs
-------
analysis/matched_pockets.csv
analysis/distance_matrix.csv
"""

from pathlib import Path
from dataclasses import dataclass
import json
import pickle

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

###############################################################################
# USER SETTINGS
###############################################################################

ALIGNMENT_FILE = Path("Analysis/01_align_structures/alignment_matrices.pkl")
POSE_JSON = Path("Data/results_KRAS/pocket_locations/pose_data.json")
POCKET_JSON = Path("Data/results_KRAS/pocket_locations/top_3_pockets.json")

OUTPUT_DIR = Path("Analysis/02_matched_pockets")
OUTPUT_DIR.mkdir(exist_ok=True)

MATCHED_OUTPUT = OUTPUT_DIR / "02_matched_pockets.csv"
DISTANCE_OUTPUT = OUTPUT_DIR / "02_distance_matrix.csv"

###############################################################################
# DATA STRUCTURE
###############################################################################

@dataclass
class Pocket:

    rank: int

    pocket_id: str

    pocket_score: float

    center: np.ndarray

    docking_score: float

###############################################################################
# LOAD INPUT FILES
###############################################################################

print("Loading alignment matrices...")

with open(ALIGNMENT_FILE, "rb") as f:
    alignments = pickle.load(f)

print("Loading top_3_pockets.json...")

with open(POCKET_JSON) as f:
    pocket_json = json.load(f)

print("Loading pose_data.json...")

with open(POSE_JSON) as f:
    pose_json = json.load(f)

###############################################################################
# HELPER FUNCTIONS
###############################################################################

def get_wildtype_key(entry):
    """
    Returns the WT key for a ligand.

    Examples
    --------
    WILD-TYE

    RAF1-WILD-TYE
    """

    for key in entry.keys():

        if "WILD" in key.upper():

            return key

    raise ValueError("Wild-type not found.")


def split_ligand_name(name):
    """
    MYC_Prednisone

    →

    MYC
    Prednisone
    """

    protein, drug = name.split("_", 1)

    return protein, drug


###############################################################################
# DOCKING PARSER
###############################################################################

def best_docking_scores(structure_dict):
    """
    Reads one structure from pose_data.json

    Parameters
    ----------
    structure_dict

    Example

    P74A_A

    Returns
    -------
    dict

    {

        1 : -58.2,

        2 : -49.1,

        3 : -63.4,

    }
    """

    scores = {}

    for pocket_number, poses in structure_dict.items():

        best = None

        for pose in poses.values():

            score = pose.get("total_score")

            if score is None:
                continue

            if best is None or score < best:

                best = score

        scores[int(pocket_number)] = best

    return scores


###############################################################################
# BUILD POCKET OBJECTS
###############################################################################

def build_pocket_list(
    pocket_list,
    docking_scores,
):
    """
    Merge

        top_3_pockets.json

    with

        pose_data.json

    into a list of Pocket objects.
    """

    pockets = []

    for pocket in pocket_list:

        rank = int(pocket["pocket_rank"])

        pockets.append(

            Pocket(

                rank=rank,

                pocket_id=pocket["selected_pocket_id"],

                pocket_score=float(

                    pocket["pocket_score"]

                ),

                center=np.array(

                    [

                        pocket["center_x"],

                        pocket["center_y"],

                        pocket["center_z"],

                    ],

                    dtype=float,

                ),

                docking_score=docking_scores.get(

                    rank,

                    np.nan,

                ),

            )

        )

    pockets.sort(

        key=lambda x: x.rank

    )

    return pockets

###############################################################################
# COORDINATE TRANSFORMATION
###############################################################################

def transform_centers(
    pockets,
    rotation,
    translation,
):
    """
    Applies the structural alignment to every
    mutant pocket center.
    """

    transformed = []

    for pocket in pockets:

        xyz = (

            pocket.center @ rotation

        ) + translation

        transformed.append(

            Pocket(

                rank=pocket.rank,

                pocket_id=pocket.pocket_id,

                pocket_score=pocket.pocket_score,

                center=xyz,

                docking_score=pocket.docking_score,

            )

        )

    return transformed

###############################################################################
# DISTANCE MATRIX
###############################################################################

def build_distance_matrix(
    wt_pockets,
    mut_pockets,
):

    wt_xyz = np.array(

        [

            p.center

            for p in wt_pockets

        ]

    )

    mut_xyz = np.array(

        [

            p.center

            for p in mut_pockets

        ]

    )

    return cdist(

        wt_xyz,

        mut_xyz,

        metric="euclidean",

    )

###############################################################################
# HUNGARIAN MATCHING
###############################################################################

def match_pockets(distance_matrix):

    rows, cols = linear_sum_assignment(

        distance_matrix

    )

    matches = []

    for r, c in zip(rows, cols):

        matches.append(

            (

                r,

                c,

                distance_matrix[r, c],

            )

        )

    return matches