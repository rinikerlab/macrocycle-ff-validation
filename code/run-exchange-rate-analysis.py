#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os.path
import json

import numpy as np
import pandas as pd

from file_structure import Simulation
from post_processing_functions import traj_position, traj_order, get_roundtrips, average_roundtrip_time, success_rate

def main():
    simulations_df = pd.read_csv("../simulations-for-exchange-rate.dat", delim_whitespace=True, names=["compound", "solvent", "method", "forcefield"])
    mkdir_if_does_not_exist("outputs")
    mkdir_if_does_not_exist("outputs/exchange-rates")
    for _, row in simulations_df.iterrows():
        sim = Simulation(row["compound"], row["forcefield"], row["method"], row["solvent"])
        print("Starting analysis on:")
        print(sim)
        run_analysis(sim)

def run_analysis(sim: Simulation):
    out = {"success_rates": success_rates(sim), "roundtrip_time": average_roundtrip_time(sim)}
    fname_base = f"outputs/exchange-rates/{sim.compound}_{sim.solvent}_{sim.method}_{sim.forcefield}"
    with open(fname_base + ".json", "w") as f:
        json.dump(out, f)

    traj_pos = traj_position(sim)
    roundtrip_events = np.array([
        get_roundtrips(traj_pos[col], sim.n_states())
        for col in traj_pos.columns
    ]).T
    np.savetxt(
        fname_base + "-roundtrip-events.csv",
        roundtrip_events,
    )

def success_rates(sim):
    order = traj_order(sim)
    rates = [success_rate(order, i) for i in range(sim.n_states()-1)]
    return rates

def mkdir_if_does_not_exist(path):
    if os.path.isdir(path):
        return
    os.mkdir(path)


if __name__ == "__main__":
    main()
