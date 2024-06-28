#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os.path

import mdtraj as md
import numpy as np
import pandas as pd
import rdkit.Chem.AllChem
from rdkit import Chem

from file_structure import Simulation
from post_processing_functions import trajectory, get_mol, extract_noe
from restraints import restraints_mol_index

def main():
    simulations_df = pd.read_csv("../simulations-only-trans.dat", delim_whitespace=True, names=["compound", "solvent", "method", "forcefield"])
    mkdir_if_does_not_exist("outputs")
    mkdir_if_does_not_exist("outputs/match-to-noes")
    for _, row in simulations_df.iterrows():
        sim = Simulation(row["compound"], row["forcefield"], row["method"], row["solvent"])
        print("Starting analysis on:")
        print(sim)
        run_analysis(sim)


def run_analysis(sim: Simulation):
    dihedrals = get_amino_acid_dihedrals(sim)
    trans = np.cos(dihedrals) < 0
    assert len(trans.shape) == 2 and trans.shape[0] >= 50_001, f"Wrong shape: {trans.shape}"
    all_trans = np.all(trans, axis=1)
    # additionally skip 10 ns as equilibration
    frame_sel = np.flatnonzero(all_trans & (np.arange(len(trans)) >= 10_000))
    df = restraints_validation_dataframe(sim, frame_sel)
    df.to_csv(f"outputs/match-to-noes/{sim.compound}_{sim.solvent}_{sim.method}_{sim.forcefield}_all-trans.csv")


def get_amino_acid_dihedrals(sim: Simulation):
    mol = get_mol(sim)
    amide_indices = get_amide_matches(mol)
    amide_index_strings = [",".join(map(str, idx)) for idx in amide_indices]

    traj = trajectory(sim)
    dihedrals = md.compute_dihedrals(traj, amide_indices)
    return dihedrals


def get_amide_matches(mol):
    query_smarts = r"[#6][NH1X3]C(=O)[#6]"
    omega_atoms = [0, 1, 2, 4]
    query = Chem.MolFromSmarts(query_smarts)
    matches = mol.GetSubstructMatches(query)
    return [tuple(match[i] for i in omega_atoms) for match in matches]


def distances_for_dataframe(traj, df):
    return [
        extract_noe(traj, row['index_a'], row['index_b'])*10.
        for _, row in df.iterrows()
    ]

def get_restraints_dataframe(simulation):
    df = restraints_mol_index(
        simulation.compound,
        simulation.solvent,
        mol=get_mol(simulation),
    )
    return df


def restraints_validation_dataframe(simulation: Simulation, frame_selection=slice(None), n_chunks=5):
    df = get_restraints_dataframe(simulation)
    traj = trajectory(simulation)[frame_selection]
    df['computed_distance'] = distances_for_dataframe(traj, df)
    # trajectory splitting
    chunks = np.array_split(np.arange(traj.n_frames), n_chunks)
    for i, indices in enumerate(chunks):
        df[f'distance_split_{i}'] = distances_for_dataframe(traj[indices], df)
    return df


def mkdir_if_does_not_exist(path):
    if os.path.isdir(path):
        return
    os.mkdir(path)


if __name__ == "__main__":
    main()
