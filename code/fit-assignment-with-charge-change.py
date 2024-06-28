#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

import rdkit.Chem.AllChem
from rdkit import Chem
import networkx

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("assignment_csv")
    parser.add_argument("neutral_mol")
    parser.add_argument("charged_mol")
    parser.add_argument("out_csv")
    args = parser.parse_args()
    neutral = Chem.MolFromMolFile(args.neutral_mol, removeHs=False)
    charged = Chem.MolFromMolFile(args.charged_mol, removeHs=False)
    print("GetSubstructMatches")
    substruct_atoms = charged.GetSubstructMatch(neutral)
    print("-- done")
    substruct_atoms = np.sort(substruct_atoms)
    # make a copy of charged without the extra hydrogens
    substruct = Chem.RWMol()
    sub_conf = Chem.Conformer()
    substruct.BeginBatchEdit()
    mapping = np.full(charged.GetNumAtoms(), -1)
    for i_atom in range(charged.GetNumAtoms()):
        if i_atom in substruct_atoms:
            substruct.AddAtom(charged.GetAtomWithIdx(i_atom))
            mapping[i_atom] = substruct.GetNumAtoms() - 1
            sub_conf.SetAtomPosition(substruct.GetNumAtoms() - 1, charged.GetConformer().GetPositions()[i_atom])
    for bond in charged.GetBonds():
        i1 = bond.GetBeginAtomIdx()
        i2 = bond.GetEndAtomIdx()
        if i1 in substruct_atoms and i2 in substruct_atoms:
            substruct.AddBond(int(mapping[i1]), int(mapping[i2]), bond.GetBondType())
    for at in substruct.GetAtoms():
        if at.GetFormalCharge() == 1:
            print("removing charge")
            at.SetFormalCharge(0)
            at.SetNumExplicitHs(0)
            at.SetNoImplicit(True)
            at.UpdatePropertyCache()
    substruct.CommitBatchEdit()
    print("AssignBondOrders")
    substruct = Chem.AllChem.AssignBondOrdersFromTemplate(neutral, substruct)
    print("-- done")
    substruct.AddConformer(sub_conf)
    Chem.AssignStereochemistryFrom3D(neutral)
    Chem.AssignStereochemistryFrom3D(substruct)
    neutral.UpdatePropertyCache()
    substruct.UpdatePropertyCache()

    # prochirality-aware structure match to find the indices in the substructure of the new mol
    sub_matching = match_atoms(substruct, neutral)
    final_matching = substruct_atoms[sub_matching]
    
    csv_atoms = pd.read_csv(args.assignment_csv, delimiter=';', names=['name', 'indices', 'comment'], dtype={'indices': str})
    csv_atoms['indices'] = [
        ','.join([str(final_matching[int(i)]) for i in row.split(',')])
        for row in csv_atoms['indices'].values
    ]
    csv_atoms.to_csv(args.out_csv, sep=';', header=False)


def match_atoms(mol_1, mol_2, conf_1=0, conf_2=0):
    """Align atom order based on the rdkit canonical rank and prochiral centers.

    When two substituents on a 4-bonded atom have the same canonical rank, uses
    a geometrical criterion to distinguish them. This way, prochiral atoms can
    be distinguished.
    This will not work with other sources of prochirality.

    Returns
    -------
    alignment: np.ndarray or None
        Array such that [mol_1.GetAtomWithIdx(i) for i in alignment] matches
        the order of mol_2.

    Raises
    ------
    ValueError: if no Isomorphism could be found.

    Examples
    --------
    >>> mol_a = Chem.AddHs(Chem.MolFromSmiles("C(C([#9])C)O"))
    >>> mol_b = Chem.AddHs(Chem.MolFromSmiles("[#9]C(C)CO"))
    >>> Chem.AllChem.EmbedMolecule(mol_a)
    0
    >>> Chem.AllChem.EmbedMolecule(mol_b)
    0
    >>> for i_a, i_b in enumerate(match_atoms(mol_a, mol_b)):
    ...     at_a = mol_a.GetAtomWithIdx(i_a)
    ...     at_b = mol_b.GetAtomWithIdx(int(i_b))
    ...     assert at_a.GetAtomicNum() == at_b.GetAtomicNum(), f"{at_a.GetAtomicNum()} != {at_b.GetAtomicNum()}"
    """
    n_atoms = mol_1.GetNumAtoms()
    assert n_atoms == mol_2.GetNumAtoms()
    graph_1 = mol_to_graph(mol_1)
    graph_2 = mol_to_graph(mol_2)
    add_tetrahedral_tiebreakers(mol_1, graph_1, conformer=conf_1)
    add_tetrahedral_tiebreakers(mol_2, graph_2, conformer=conf_2)
    iso = networkx.vf2pp_isomorphism(graph_1, graph_2, node_label="label")
    if iso is None:
        raise ValueError("Graph isomerism between molecules failed. Are the bonds equivalent?")
    return np.array([iso[i] for i in range(n_atoms)])


def group(values):
    "Return a dict of value -> list of indices"
    out = {}
    for i, val in enumerate(values):
        out.setdefault(val, []).append(i)
    return out


def get_plane_normal(x1, x2, x3):
    "Return a unit array normal to the plane defined by x1, x2, x3"
    x1 = np.array(x1)
    x2 = np.array(x2)
    x3 = np.array(x3)
    normal = np.cross(x1-x2, x3-x2)
    return normal / np.linalg.norm(normal)


def mol_to_graph(mol: Chem.Mol):
    """Networkx graph where the rdkit canonical rank is used as Node label."""
    graph = networkx.Graph()
    canon_rank = Chem.CanonicalRankAtoms(mol, breakTies=False)
    for i, atom in enumerate(mol.GetAtoms()):
        graph.add_node(i, label=(canon_rank[i],), elem=atom.GetAtomicNum())
    for i, bond in enumerate(mol.GetBonds()):
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return graph


def add_tetrahedral_tiebreakers(mol, graph, conformer=0):
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    conf = mol.GetConformer(conformer)
    xyz = np.array(conf.GetPositions())
    for i_atom in range(mol.GetNumAtoms()):
        neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(i_atom).GetNeighbors()]
        nbr_ranks = group([ranks[i] for i in neighbors])
        # 4-bonded atom where two have the same rank
        if len(neighbors) == 4 and len(nbr_ranks) == 3:
            prochiral_nbrs = next(group for group in nbr_ranks.values() if len(group) == 2)
            prochiral = [neighbors[i] for i in prochiral_nbrs]
            one_atom_groups = sorted([i for i, group in nbr_ranks.items() if len(group) == 1])
            unique_nbrs = [nbr_ranks[i][0] for i in one_atom_groups]
            unique = [neighbors[i] for i in unique_nbrs]
            plane_normal = get_plane_normal(xyz[unique[0]], xyz[i_atom], xyz[unique[1]])
            for i_pro in prochiral:
                direction = xyz[i_pro] - xyz[i_atom]
                prochirality = int(np.sum(direction * plane_normal) > 0)
                graph.nodes[i_pro]['label'] = graph.nodes[i_pro]['label'] + (prochirality,)

if __name__ == "__main__":
    main()
