#! /bin/bash

. ./profile.inc

for YEAR in `seq 2019 2022`; do

CHECKDIR=/g100_scratch/userexternal/gbolzon0/V9C/transition/wrkdir/DA/${YEAR}
DA_DIR=/g100_scratch/userexternal/gbolzon0/V9C/transition/wrkdir/MODEL/DA__FREQ_1/${YEAR}
CHECKDATA=/g100_scratch/userexternal/gbolzon0/V9C/transition/wrkdir/POSTPROC/validation/post_check_data
IMG_DIR=/g100_scratch/userexternal/gbolzon0/V9C/transition/wrkdir/POSTPROC/validation/IMG/post_check

mkdir -p $CHECKDATA $IMG_DIR

#opa_prex_or_die "python post_check_float.py -i $CHECKDIR -d  $DA_DIR -o $CHECKDATA "

done



opa_prex_or_die "python plot_post_check_float.py -i $CHECKDATA -o $IMG_DIR"