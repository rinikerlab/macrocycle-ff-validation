#!/bin/bash

cd outputs || exit
acpype -p prepared.parm7 -x prepared.rst7 -b converted
