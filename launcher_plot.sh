#! /bin/bash

source $HOME/sequence.sh

export PYTHONPATH=$PWD:$PYTHONPATH
python matchup_seasonal_WE.py  -v Ed380f
python matchup_seasonal_WE.py  -v Ed412f
python matchup_seasonal_WE.py  -v Ed490f
