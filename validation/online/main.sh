#! /bin/bash

#OPA_RUNDATE has to be known at this level

#step 0  -- set of starttime, endtime
STARTTIME_a=20160305
END__TIME_a=20160308

STARTTIME_f=$END__TIME_a
END__TIME_f=20160312

# source bit.sea/config.sh is needed 

# step 1 -- unzip from archive
ARCHIVE_DIR=/pico/home/usera07ogs/a07ogs00/OPA/V4/archive # deve diventare V2

ONLINE_VALIDATION_DIR=/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data_02
TMP_DIR=${ONLINE_VALIDATION_DIR}/TMP
PROFILERDIR=${ONLINE_VALIDATION_DIR}/PROFILATORE
IMAGES_DIR=./outdir/
MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc

mkdir -p $ONLINE_VALIDATION_DIR
python archive_extractor.py --type analysis -st ${STARTTIME_a} -et ${END__TIME_a}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}
python archive_extractor.py --type forecast -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}


#step 2 -- re-do Forcings generation
./forcings_gen.sh -d ${ONLINE_VALIDATION_DIR}

# step 3 -- aggregate variables to have chl and all variables (physical and biological) on the same ave files

python var_aggregator.py --biodir ${ONLINE_VALIDATION_DIR}/output_bio --physdir ${ONLINE_VALIDATION_DIR}/FORCINGS -t $TMP_DIR

# step 4
python float_extractor.py -st ${STARTTIME} -et ${END__TIME} -i $TMP_DIR -b $PROFILERDIR  -o $IMAGES_DIR -m $MASKFILE
