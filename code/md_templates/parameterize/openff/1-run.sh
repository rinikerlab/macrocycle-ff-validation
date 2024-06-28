#!/bin/bash

cp ../../create-initial-structure/outputs/rdk_model.mol . || exit
mkdir -p outputs || exit
( cd outputs; obabel -imol ../rdk_model.mol -xn -opdb -O model.pdb )
python parameterize.py
