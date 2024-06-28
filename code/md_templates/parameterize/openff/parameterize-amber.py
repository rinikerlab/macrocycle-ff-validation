#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

from openff.toolkit import Molecule
from openff.toolkit import ForceField
from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--molfile', default='rdk_model.mol')
    parser.add_argument('--top_out', default='outputs-amber/model.parm7')
    parser.add_argument('--crd_out', default='outputs-amber/model.inpcrd')
    parser.add_argument('--infile_out', default='outputs-amber/example.in')
    if not os.path.exists('outputs-amber'):
        os.mkdir('outputs-amber')
    args = parser.parse_args()
    mol = Molecule.from_file(args.molfile)
    topology = mol.to_topology()
    forcefield = ForceField("openff-2.0.0.offxml")
    interchange = Interchange.from_smirnoff(
        force_field=forcefield,
        topology=topology,
    )
    interchange.to_prmtop(args.top_out)
    interchange.to_inpcrd(args.crd_out)
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.periodic = True
    mdconfig.write_sander_input_file(input_file=args.infile_out)

if __name__ == "__main__":
    main()
