#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path

import rdkit.Chem.AllChem
from rdkit import Chem

import numpy as np
import mdtraj as md
import pandas as pd

from file_structure import Simulation
from post_processing_functions import trajectory, get_mol

def main():
    simulations_df = pd.read_csv("../simulations-for-cis-fraction.dat", delim_whitespace=True, names=["compound", "solvent", "method", "forcefield"])
    mkdir_if_does_not_exist("outputs")
    mkdir_if_does_not_exist("outputs/cis-fraction")
    for _, row in simulations_df.iterrows():
        sim = Simulation(row["compound"], row["forcefield"], row["method"], row["solvent"])
        print("Starting analysis on:")
        print(sim)
        run_analysis(sim)


def run_analysis(sim: Simulation):
    mol = get_mol(sim)
    amide_indices = get_amide_matches(mol)
    amide_index_strings = [",".join(map(str, idx)) for idx in amide_indices]

    traj = trajectory(sim)
    dihedrals = md.compute_dihedrals(traj, amide_indices)
    df = pd.DataFrame(dihedrals, columns=amide_index_strings)
    fname_base = f"outputs/cis-fraction/{sim.compound}_{sim.solvent}_{sim.method}_{sim.forcefield}"
    df.to_csv(fname_base + "_full.csv")
    cis_fraction = (np.cos(df) > 0).mean().rename("cis_fraction").rename_axis("dihedral")
    cis_fraction.to_csv(fname_base + "_summary.csv")


def get_amide_matches(mol):
    query_smarts = r"[#6][NH1X3]C(=O)[#6]"
    omega_atoms = [0, 1, 2, 4]
    query = Chem.MolFromSmarts(query_smarts)
    matches = mol.GetSubstructMatches(query)
    return [tuple(match[i] for i in omega_atoms) for match in matches]


def mkdir_if_does_not_exist(path):
    if os.path.isdir(path):
        return
    os.mkdir(path)


if __name__ == "__main__":
    main()
