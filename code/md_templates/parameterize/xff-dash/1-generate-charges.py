#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from rdkit import Chem
from serenityff.charge.tree.dash_tree import DASHTree

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles', default='../../smiles.smi')
    parser.add_argument('--out', default='charges.dat')
    args = parser.parse_args()
    with open(args.smiles) as f:
        smiles = next(f).strip()
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    tree = DASHTree(num_processes=1)
    tree_output = tree.get_molecules_partial_charges(mol)
    charges = tree_output['charges']
    max_stdev = max(tree_output['std'])
    print("Max standard deviation of charges:", max_stdev)
    with open(args.out, "w") as f:
        print(*charges, sep="\n", file=f)


if __name__ == "__main__":
    main()
