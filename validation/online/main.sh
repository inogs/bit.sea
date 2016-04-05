#! /bin/bash

#OPA_RUNDATE has to be known at this level

#step 0  -- set of starttime, endtime
STARTTIME=20160301
END__TIME=20160308

# source bit.sea/config.sh is needed 

# step 1 -- unzip from archive
ARCHIVE_DIR=/pico/home/usera07ogs/a07ogs00/OPA/V4/archive # deve diventare V2

ONLINE_VALIDATION_DIR=/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data
UNZIPPED_AVE_DIR=${ONLINE_VALIDATION_DIR}/
TMP_DIR=${ONLINE_VALIDATION_DIR}/TMP
PROFILERDIR=${ONLINE_VALIDATION_DIR}/PROFILATORE
IMAGES_DIR=./outdir/
MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc


python archive_extractor.py  -st ${STARTTIME} -et ${END__TIME} -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}

#step 2 -- re-do Forcings generation
./forcings_gen.sh -d ${ONLINE_VALIDATION_DIR}


# step 3 -- aggregate variables to have chl and all variables (physical and biological) on the same ave files

python var_aggregator.py --biodir ${ONLINE_VALIDATION_DIR}/output_bio --physdir ${ONLINE_VALIDATION_DIR}/FORCINGS -t $TMP_DIR

# step 4
python float_extractor.py -st ${STARTTIME} -et ${END__TIME} -i $TMP_DIR -b $PROFILERDIR  -o $IMAGES_DIR -m $MASKFILE
