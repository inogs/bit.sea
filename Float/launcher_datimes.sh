#! /bin/bash

export ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/ONLINE_V5C/

DAFREQ=1
DA_HH=13

DATIMES_DIR=./
OUTFILE=${DATIMES_DIR}/daTimes

for VARDA in N3n P_l; do
python profile_dates_DAfreq.py -f $DAFREQ -t $DA_HH -s 20170102 -e 20170201 -v ${VARDA} -o ${DATIMES_DIR}/daTimes_float_${VARDA}
done

# Modify varlist in merge_daTimes.py
DA_TIMES_SAT=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_04/wrkdir/MODEL/daTimes_sat
python merge_daTimes.py -s $DA_TIMES_SAT -d $DATIMES_DIR -o ${OUTFILE}
