#!/usr/bin/env python3
"""
alignment.py

Utility functions for structural alignment of AlphaFold models.

Author: Leo + ChatGPT

Requirements
------------
biopython
numpy
"""

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from Bio.PDB import PDBParser
from Bio.PDB.Superimposer import Superimposer


###############################################################################
# PDB LOADING
###############################################################################

parser = PDBParser(QUIET=True)


def load_structure(pdb_file: str):
    """
    Load a PDB structure.

    Parameters
    ----------
    pdb_file : str

    Returns
    -------
    Bio.PDB.Structure.Structure
    """
    return parser.get_structure(Path(pdb_file).stem, pdb_file)


###############################################################################
# CA EXTRACTION
###############################################################################

def extract_ca_atoms(structure):
    """
    Extract all CA atoms from chain A.

    Returns
    -------
    dict

    {
        residue_number : CA_atom
    }
    """

    model = next(structure.get_models())

    if "A" not in model:
        raise ValueError("Chain A not found.")

    chain = model["A"]

    atoms = {}

    for residue in chain:

        if residue.id[0] != " ":
            continue

        if "CA" not in residue:
            continue

        resnum = residue.id[1]

        atoms[resnum] = residue["CA"]

    return atoms


###############################################################################
# MUTATION PARSING
###############################################################################

def parse_mutation(protein: str, folder_name: str):
    """
    Converts model folder names to the canonical mutation name.

    Examples
    --------
    protein="MYC", folder="pP74A"
        -> MYC-P74A

    protein="ACVR1", folder="pG328E"
        -> ACVR1-G328E

    protein="BRAF", folder="BRAF-V600E"
        -> BRAF-V600E

    protein="RAF1", folder="RAF1-S257L"
        -> RAF1-S257L

    wild-type
        -> None
    """

    name = folder_name

    if "WILD" in name.upper():
        return None

    if name.lower() == "wild-type":
        return None

    if name.startswith("p"):
        name = name[1:]

    if not name.startswith(f"{protein}-"):
        name = f"{protein}-{name}"

    return name


###############################################################################
# MATCH RESIDUES
###############################################################################

def match_ca_atoms(
    wt_atoms: Dict,
    mut_atoms: Dict,
    mutation: str = None,
):
    """
    Match CA atoms using residue numbering.

    Parameters
    ----------
    wt_atoms

    mut_atoms

    mutation
        Example: P74A

    Returns
    -------
    wt_list
    mut_list
    aligned_residue_numbers
    """

    common = sorted(
        set(wt_atoms.keys()) &
        set(mut_atoms.keys())
    )

    if mutation is not None:

        number = ""

        for c in mutation:

            if c.isdigit():
                number += c

        if number != "":

            mutation_position = int(number)

            common = [
                x
                for x in common
                if x != mutation_position
            ]

    wt_list = []
    mut_list = []

    for residue in common:

        wt_list.append(wt_atoms[residue])

        mut_list.append(mut_atoms[residue])

    return wt_list, mut_list, common


###############################################################################
# SUPERPOSITION
###############################################################################

def superimpose_structures(
    wt_pdb: str,
    mut_pdb: str,
    mutation: str = None,
):
    """
    Superimpose mutant onto WT.

    Returns
    -------
    dict
    """

    wt_structure = load_structure(wt_pdb)

    mut_structure = load_structure(mut_pdb)

    wt_atoms = extract_ca_atoms(wt_structure)

    mut_atoms = extract_ca_atoms(mut_structure)

    wt_list, mut_list, residues = match_ca_atoms(
        wt_atoms,
        mut_atoms,
        mutation,
    )

    sup = Superimposer()

    sup.set_atoms(
        wt_list,
        mut_list,
    )

    rotation, translation = sup.rotran

    rmsd = sup.rms

    return {

        "rotation": rotation,

        "translation": translation,

        "rmsd": rmsd,

        "aligned_residues": len(residues),

        "residue_numbers": residues,

        "wt_atoms": wt_list,

        "mut_atoms": mut_list,
    }


###############################################################################
# COORDINATE TRANSFORMATION
###############################################################################

def transform_coordinates(
    coordinates: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
):
    """
    Transform coordinates from mutant reference frame
    into WT reference frame.

    Parameters
    ----------
    coordinates : (N,3)

    Returns
    -------
    transformed : (N,3)
    """

    coordinates = np.asarray(coordinates)

    return np.dot(
        coordinates,
        rotation
    ) + translation