#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import rdkit.Chem.AllChem
from rdkit import Chem

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_file", default="../smiles.smi", required=False)
    parser.add_argument("--out", default="outputs/rdk_model.mol", required=False)
    args = parser.parse_args()
    with open(args.smiles_file) as f:
        smiles = next(f).strip()
    mol = Chem.MolFromSmiles(smiles)
    with_hydrogen = Chem.AddHs(mol)
    params = Chem.AllChem.ETKDGv3()
    params.randomSeed = 0xf00d
    conf_id = Chem.AllChem.EmbedMolecule(with_hydrogen, params=params)
    assert conf_id == 0
    # order = get_atom_order(with_hydrogen)
    # renumbered = Chem.RenumberAtoms(with_hydrogen, order)
    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    Chem.MolToMolFile(with_hydrogen, args.out)

def get_atom_order(mol):
    bonds = {i: [] for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        bonds[bond.GetBeginAtomIdx()] = bond.GetEndAtomIdx()
        bonds[bond.GetEndAtomIdx()] = bond.GetBeginAtomIdx()
    def atom_sort_key(index):
        atom = mol.GetAtomWithIdx(index)
        if atom.GetAtomicNum() == 1:
            return (bonds[atom.GetIdx()], atom.GetIdx())
        return (atom.GetIdx(),)
    return sorted(range(mol.GetNumAtoms()), key=atom_sort_key)

if __name__ == "__main__":
    main()
