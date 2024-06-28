import rdkit.Chem.AllChem
from rdkit import Chem
import numpy as np
import networkx


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
    assert n_atoms == mol_2.GetNumAtoms(), f"Different number of atoms: {n_atoms} != {mol_2.GetNumAtoms()}"
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


def test_match_atoms():
    mol1 = Chem.AddHs(Chem.MolFromSmiles('C1C(C)=CC=C1'))
    Chem.AllChem.EmbedMolecule(mol1)
    conf1 = mol1.GetConformer(0)
    mol2 = Chem.RWMol()
    # mol2 will be mol1, but with the 2 prochiral hydrogens replaced.
    replacement = {6: 7, 7: 6}
    atoms = mol1.GetAtoms()
    n_atoms = len(atoms)
    conf2 = Chem.Conformer()
    for i in range(n_atoms):
        mol2.AddAtom(atoms[replacement.get(i, i)])
        conf2.SetAtomPosition(i, conf1.GetAtomPosition(replacement.get(i, i)))
    for bond in mol1.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        mol2.AddBond(replacement.get(begin, begin), replacement.get(end, end))
    conf_id = mol2.AddConformer(conf2)
    match = match_atoms(mol1, mol2, conf_2=conf_id)
    assert list(match) == [0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 13]
