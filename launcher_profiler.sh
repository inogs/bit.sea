#! /bin/bash

source $HOME/sequence.sh

export MASKFILE=/galileo/home/userexternal/eterzic0/MASKFILE_70/meshmask.nc


export PYTHONPATH=$PYTHONPATH:$PWD

python profiler_ET.py


