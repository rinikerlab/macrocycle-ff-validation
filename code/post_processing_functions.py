from functools import cache
import re
import gzip
import warnings

import numpy as np
import mdtraj as md
import pandas as pd
import rdkit.Chem.AllChem
from rdkit import Chem
from file_structure import Simulation, molfile_path


# Load

def load_pdb_and_assign_bonds(fname, refmol) -> Chem.Mol:
    mol = Chem.MolFromPDBFile(fname, removeHs=False)
    with_bonds = Chem.AllChem.AssignBondOrdersFromTemplate(refmol, mol)
    Chem.AllChem.AssignStereochemistryFrom3D(with_bonds)
    return with_bonds


def get_mol(simulation) -> Chem.Mol:
    refmol = Chem.MolFromMolFile(molfile_path(simulation.compound), removeHs=False)
    fname = simulation.pdb_filename()
    return load_pdb_and_assign_bonds(fname, refmol)


def load(path, top, make_whole=True) -> md.Trajectory:
    "Load a trajectory into MDTraj and assign bonds"
    traj: md.Trajectory = md.load(path, top=top)
    if make_whole:
        if traj.unitcell_lengths is not None:
            # if len(traj.top._bonds) == 0:
            #     bonds = guess_bonds(traj)
            #     traj.make_molecules_whole(sorted_bonds=bonds)
            # else:
                traj.make_molecules_whole()
        else:
            warnings.warn("Cannot make molecules whole because the trajectory has no unitcell")
    return traj


def trajectory(simulation: Simulation) -> md.Trajectory:
    return load(simulation.trajectory_path(), topology(simulation))


def load_pdbfile(fname) -> md.Trajectory:
    return md.load(fname)


def pdbfile(simulation: Simulation) -> md.Trajectory:
    return load_pdbfile(simulation.pdb_filename())


def topology(simulation: Simulation) -> md.Topology:
    return pdbfile(simulation).top


# Parse Gromacs logfiles with plumed exchange info

@cache
def exchanges_re(n_states):
    """
    Examples
    --------
    >>> test_str = '0    1    2 x  3    4 x  5    6 x  7 x  8    9 x 10   11'
    >>> print(exchanges_re(12).findall(test_str))
    [('', '', 'x', '', 'x', '', 'x', 'x', '', 'x', '')]
    """
    return re.compile('  *(x?)  *'.join([str(i) for i in range(n_states)]))


def parse_exchanges(line, n_states):
    """
    Examples
    --------
    >>> test_str = '0    1    2 x  3    4 x  5    6 x  7 x  8    9 x 10   11'
    >>> print(parse_exchanges(test_str, 12))
    [(2, 3), (4, 5), (6, 7), (7, 8), (9, 10)]
    """
    x_or_empty = exchanges_re(n_states).findall(line.strip())[0]
    assert len(x_or_empty) == n_states - 1, f"wrong number of matches: {len(x_or_empty)}"
    exchanges = []
    for i in range(n_states - 1):
        if x_or_empty[i] == 'x':
            exchanges.append((i, i+1))
    return exchanges


def read_trajectory_timeline(fh, n_states):
    """Parse the traj_order from a Gromacs+plumed logfile.
    See traj_order for an explation of the output data"""
    trajectories = list(range(n_states))
    out = [list(trajectories), list(trajectories)]
    for line in fh:
        if line.startswith('Repl ex '):
            exchange_string = line[len('Repl ex '):].strip()
            for s1, s2 in parse_exchanges(exchange_string, n_states):
                trajectories[s1], trajectories[s2] = trajectories[s2], trajectories[s1]
            out.append(list(trajectories))
    return out


# def read_exchange_times(fh):
#     out = [0., 0.]  # prepend 0 for consistency with read_trajectory_timeline
#     for line in fh:
#         if line.startswith('Replica exchange at step '):
#             out.append(format(float(line.split()[-1]), '.1f'))
#     return out


def read_traj_order(fh, n_states):
    return pd.DataFrame(read_trajectory_timeline(fh, n_states)).rename_axis('frame').rename_axis('position', axis=1)


def traj_order(simulation: Simulation) -> pd.DataFrame:
    """Which continious trajectory was in the Nth ensemble 
    at each frame?
    
    Returns a Dataframe with shape n_frames, n_simulations.
    
    E.g., if the 2nd continuous trajectory (index 1) was
    in the sixth ensemble (index 5) at frame 1000, then
    traj_order(sim).loc[1000, 5] == 1
    """
    with gzip.open(simulation.logfile_path(), 'rt') as fh:
        return read_traj_order(fh, simulation.n_states())

    
# def read_times(fh):
#     return pd.Series(read_exchange_times(fh)).rename_axis('frame').rename('time')


# @lru_cache(2)
# def times(simulation: Simulation):
#     with open(simulation.logfile_path()) as fh:
#         return read_times(fh)


def traj_position(simulation: Simulation) -> pd.DataFrame:
    """In which position (ensemble) was the Nth continuous 
    trajectory at each frame?
    
    Returns a Dataframe with shape n_frames, n_simulations.
    
    E.g., if the 2nd continuous trajectory (index 1) was
    in the sixth ensemble (index 5) at frame 1000, then
    traj_position(sim).loc[1000, 1] == 5
    """
    order = traj_order(simulation)
    return pd.DataFrame(np.argsort(order), index=order.index).rename_axis('trajectory', axis=1)


# Demux

def demux(trajectories, trajectory_positions: pd.DataFrame):
    """Demux replica exchange trajectories
    
    Given a set of trajectories that are continuous in the replica,
    but non-continuous in the coordinate space, creates a set of
    trajectories that are continous in coordinate space instead.
    
    Those trajectories can be used e.g., to unwrap the trajectory,
    since they don't contain large jumps, but they should not be
    used for analysis since they are not drawn from a specific
    ensemble.
    
    trajectory_positions should be the dataframe output by
    traj_position.
    """
    n_frames = min(t.n_frames for t in trajectories)
    if n_frames > len(trajectory_positions):
        print(f"WARNING: shortest trajectory has {n_frames} frames, which is more than len(trajectory_positions) ({len(trajectory_positions)}). Truncating!")
        #trajectory_positions = trajectory_positions.iloc[:n_frames]
        n_frames = len(trajectory_positions)
    assert n_frames <= len(trajectory_positions), f"{n_frames} > {len(trajectory_positions)}"
    frames = np.arange(n_frames)
    out = [t[:n_frames] for t in trajectories]
    all_xyz = np.stack([t.xyz for t in out])
    for i in range(len(trajectories)):
        out[i].xyz = all_xyz[trajectory_positions[i].iloc[:n_frames].values.ravel(), frames, :]
    return out


def demuxed_trajectories(simulation: Simulation, stride: int):
    sims = [simulation.with_param("state", i) for i in range(simulation.n_states())]
    return demux([trajectory(s) for s in sims], traj_position(simulation)[::stride])


# Roundtrip-time related

def get_roundtrips(arr, n_states):
    """
    Given an array with the state (ensemble) in every frame of a
    continuous trajectory from replica exchange, return an array
    of roundtrip event types
    (0: no roundtrip; -1: reached the bottom; 1: reached the top)
    
    Examples
    --------
    >>> get_roundtrips([0, 1, 2, 1, 2, 0], 3)
    array([ 0,  0,  1,  0,  0, -1])
    """
    assert len(np.shape(arr)) == 1
    direction = "unknown"
    roundtrip_type = np.zeros(len(arr), dtype=int)
    for i, state in enumerate(arr):
        if state == 0:
            if direction == "down":
                roundtrip_type[i] = -1
            direction = "up"
        elif state == n_states - 1:
            if direction == "up":
                roundtrip_type[i] = 1
            direction = "down"
    return roundtrip_type


def count_roundtrips(arr, n_states):
    return np.sum(get_roundtrips(arr, n_states=n_states) != 0) / 2


def average_roundtrip_time(sim: Simulation):
    traj = trajectory(sim)
    # convert from ps to ns
    simulation_time = (traj.time[-1] - traj.time[0]) / 1000
    print(f"{simulation_time=}")
    n_roundtrips = traj_position(sim).apply(lambda x: count_roundtrips(x, sim.n_states()))
    roundtrip_time = simulation_time / np.mean(n_roundtrips)
    return roundtrip_time


# Compute the exchange rate between replica exchange ensembles

def get_attempts(traj_order, n_state):
    """Check the success of exchanges between n_state and n_state+1.
    
    Returns a binary array of length n_frames//2.
    Each value represents an exchange attempt between n_state and n_state + 1, and is True if it succeeded.
    """
    state_arr = traj_order[n_state].values.ravel()
    # First attempts to switch even states "up" and uneven ones "down", then the other way around
    # skip first frame for even states
    state_arr = state_arr[(n_state+1) % 2:]
    n_attempts = len(state_arr) // 2
    attempts = state_arr[:n_attempts*2].reshape(-1, 2)
    #print(attempts[:10])
    return attempts[:, 0] != attempts[:, 1]


def success_rate(traj_order, n_state):
    attempts = get_attempts(traj_order, n_state)
    return np.sum(attempts) / len(attempts)


# NOE-related

def extract_noe(traj: md.Trajectory, atoms1, atoms2, exponent=-6):
    if isinstance(atoms1, int):
        atoms1 = [atoms1]
    if isinstance(atoms2, int):
        atoms2 = [atoms2]
    x1 = traj.atom_slice(atoms1).xyz
    x2 = traj.atom_slice(atoms2).xyz
    dx = x1[:, :, np.newaxis, :] - x2[:, np.newaxis, :, :]
    distance = np.linalg.norm(dx, axis=-1)
    return noe_avg(distance, exponent=exponent)

def noe_avg(a, exponent=-6, axis=None):
    a = np.asarray(a, dtype=float)
    return np.mean(a**exponent, axis=axis)**(1/exponent)

assert noe_avg([1, 2]) == np.mean([1, 2**(-6)])**(-1/6), noe_avg([1, 2])
