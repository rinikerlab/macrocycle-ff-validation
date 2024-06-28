#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from openff.toolkit import Molecule
from openff.toolkit import ForceField
from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--molfile', default='rdk_model.mol')
    parser.add_argument('--top_out', default='outputs/model.top')
    parser.add_argument('--gro_out', default='outputs/model.gro')
    parser.add_argument('--mdp_out', default='outputs/example.mdp')
    # parser.add_argument('--smiles', default=r'COC(=O)[C@H]1CSCC2=C(C=CC=C2)C(=O)OC[C@H](NC(=O)[C@@H]2CCCN2C(C)=O)C(=O)N1')
    args = parser.parse_args()
    mol = Molecule.from_file(args.molfile)
    # pdbfile = app.PDBFile(args.pdbfile)
    topology = mol.to_topology()
    # off_topology = Topology.from_openmm(
        # omm_topology, unique_molecules=[mol]
    # )
    forcefield = ForceField("openff-2.0.0.offxml")
    interchange = Interchange.from_smirnoff(
        force_field=forcefield,
        topology=topology,
    )
    interchange.to_top(args.top_out)
    interchange.to_gro(args.gro_out)
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.periodic = True
    mdconfig.write_mdp_file(mdp_file=args.mdp_out)

if __name__ == "__main__":
    main()
