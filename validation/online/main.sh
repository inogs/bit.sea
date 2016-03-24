#! /bin/bash

#OPA_RUNDATE has to be known at this level

# step 1 -- unzip from archive
python archive_extractor.py 

#step 2 -- re-do Forcings generation

./forcings_gen.sh 


# step 3 -- aggregate variables to have chl


# step 4
python  float_extractor.py 
