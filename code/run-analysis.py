#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run the analysis to get violation plot data from a REST2 simulation. By default
(if not parameters are given), this expects a file called "../simulations.dat",
and gets the list of simulation folders from there. You can also specify a
single simulation using the commandline arguments.
"""

import os.path
import argparse

import numpy as np
import pandas as pd
from file_structure import Simulation
from post_processing_functions import trajectory, get_mol, extract_noe
from restraints import restraints_mol_index

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compound")
    parser.add_argument("--forcefield")
    parser.add_argument("--method")
    parser.add_argument("--solvent")
    args = parser.parse_args()

    # if called without any arguments
    if args.compound is None and args.forcefield is None and args.method is None and args.solvent is None:
        run_on_all_default()
    else:
        sim = Simulation(args.compound, args.forcefield, args.method, args.solvent)
        run_analysis(sim)


def run_on_all_default():
    """Run on all the simulations listed in simulations.dat (csv file)"""
    simulations_df = pd.read_csv("../simulations.dat", delim_whitespace=True, names=["compound", "solvent", "method", "forcefield"])
    mkdir_if_does_not_exist("outputs")
    mkdir_if_does_not_exist("outputs/match-to-noes")
    for _, row in simulations_df.iterrows():
        sim = Simulation(row["compound"], row["forcefield"], row["method"], row["solvent"])
        print("Starting analysis on:")
        print(sim)
        run_analysis(sim)


def run_analysis(sim: Simulation):
    frame_sel = slice(10_000, 50_000)
    df = restraints_validation_dataframe(sim, frame_sel)
    df.to_csv(f"outputs/match-to-noes/{sim.compound}_{sim.solvent}_{sim.method}_{sim.forcefield}.csv")


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
