#! /bin/bash

#OPA_RUNDATE has to be known at this level

# step 1 -- unzip from archive
python archive_extractor.py  -st 20160301 -et 20160308 -adir /pico/home/usera07ogs/a07ogs00/OPA/V4/archive/ -odir /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/

#step 2 -- re-do Forcings generation
./forcings_gen.sh -d /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data


# step 3 -- aggregate variables to have chl
python var_aggregator.py --biodir /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/output_bio --physdir /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/FORCINGS/ -t /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/TMP

# step 4
python float_extractor.py -st 20160301 -et 20160308 -idir /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/TMP/ -bdir /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/PROFILATORE/ -odir ./outdir/
