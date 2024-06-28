import pandas as pd
from rdkit import Chem

from atom_matching import match_atoms
from file_structure import molfile_path, restraintfile_path, atom_assignment_path


def restraints(compound: str, solvent: str) -> pd.DataFrame:
    return (
        pd.read_csv(restraintfile_path(compound, solvent), dtype={"proton_a": str, "proton_b": str})
        .astype({'proton_a': str, 'proton_b': str, 'restraint_number': str})
    )


def restraints_mol_index(compound: str, solvent: str, mol=None, use_ring_distance=False, drop_ref=True) -> pd.DataFrame:
    ref_mol = Chem.MolFromMolFile(molfile_path(compound), removeHs=False)
    ref_assignment = h_assignment(compound, solvent)
    restraints_df = restraints(compound, solvent)
    return restraints_mol_index_df(restraints_df, ref_assignment, ref_mol, mol=mol, use_ring_distance=use_ring_distance, drop_ref=drop_ref)


def restraints_mol_index_df(restraints: pd.DataFrame, ref_assignment: pd.Series, ref_mol: Chem.Mol, mol=None, use_ring_distance=False, drop_ref=True) -> pd.DataFrame:
    mol = mol or ref_mol
    alignment = match_atoms(ref_mol, mol)
    assign = pd.Series(
        [[alignment[i] for i in indices] for indices in ref_assignment.values],
        index=ref_assignment.index,
    )
    df = (
        restraints
        .join(assign.rename('index_a'), on='proton_a')
        .join(assign.rename('index_b'), on='proton_b')
    )
    nan_df = df.loc[pd.isna(df['index_a']) | pd.isna(df['index_b'])]
    if len(nan_df) > 0:
        print("Warning: found NaN entries in the DataFrame:")
        print(nan_df)
    atoms1 = [atoms[0] for atoms in df['index_a'].values]
    atoms2 = [atoms[0] for atoms in df['index_b'].values]
    df['top_distance'] = [get_bonds_distance(mol, a1, a2, in_ring=use_ring_distance) for a1, a2 in zip(atoms1, atoms2)]
    if drop_ref:
        df = drop_reference_rows(df)
    return df


def drop_reference_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Remove restraints named Ref., or between an atom with itself.

    It can happen that there is a restraint with an atom and itself, if the
    prochiral assignment is incomplete and the restraint is between hydrogen
    atoms attached to the same heavy atom.
    """
    df = df.loc[~df['restraint_number'].str.startswith("Ref")]
    df = df.query('index_a != index_b')
    return df


def get_bonds_distance(mol, atom1, atom2, in_ring=False):
    if in_ring:
        return get_bonds_distance_in_ring(mol, atom1, atom2)
    return Chem.GetDistanceMatrix(mol)[atom1, atom2]


def get_bonds_distance_in_ring(mol, atom1, atom2):
    if atom1 == atom2:
        return 0
    ring = set(max(Chem.GetSSSR(mol), key=len))
    path = Chem.GetShortestPath(mol, atom1, atom2)
    return len([a for a in path if a in ring])


def h_assignment(compound, solvent=None):
    """
    Get a hydrogen assignment and put it into a standard format.

    The format is:
    * type: pd.Series
    * data is lists of h_indices for each named atom in the NOE list
    * the keys are the atoms names as they occur in the NOE list
    * the name of the series contains metadata: compound_solvent.
      Note that the solvent might be the default-solvent, if no specific
      versions for the given inputs exist.
    """
    return parse_h_assignment(atom_assignment_path(compound, solvent))


def parse_h_assignment(file):
    df = pd.read_csv(
        file,
        names=['h_name', 'atom_index', 'notes'],
        delimiter=";",
        dtype=str,
    ).set_index("h_name")  # don't use index_col, otherwise the index might be converted to int.
    h_assign = df['atom_index'].map(lambda x: [int(elem) for elem in x.split(",")])
    return h_assign


def test_h_assignment():
    from io import StringIO
    test_file = StringIO("""H1;1,2;
H2;2,3;test""")
    output = parse_h_assignment(test_file)
    assert all(output.index == ["H1", "H2"])
    assert output["H1"] == [1, 2]
    assert output["H2"] == [2, 3]

def test_h_assignment_single_col():
    from io import StringIO
    test_file = StringIO("""H1;1;
H2;2;test""")
    output = parse_h_assignment(test_file)
    assert all(output.index == ["H1", "H2"])
    assert output["H1"] == [1]
    assert output["H2"] == [2]
    assert len(output) == 2
