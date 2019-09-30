#! /bin/bash

STARTTIME_f=20190820
END__TIME_f=20190915
TMP_DIR__ACT=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/AVE/ANALYSIS
PROFILERDIRA=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/PROFILATORE_7VARS/
IMG_DIR__ACT=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/
export MASKFILE=/gpfs/work/OGS_prod_0/OPA/V5C/prod/wrkdir/2/MODEL/meshmask.nc
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i $TMP_DIR__ACT -b $PROFILERDIRA  -o $IMG_DIR__ACT -m $MASKFILE
