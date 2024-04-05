#! /bin/bash

source ../../Sat/profile.inc
source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence3.sh


TRAIN_DIR=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/Float/ppcon/results #inputs



export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C/
export ONLINE_REPO=/g100_scratch/userexternal/gbolzon0/carolina/ONLINE_REPO/
ONLINE_REPO_CLUSTERING=${ONLINE_REPO}/PPCON/clustering/

mkdir -p $ONLINE_REPO_CLUSTERING

FILE_INPUT=/g100_work/OGS_prodC/OPA/V10C-prod/log/daily/DIFF_floats.20240404-22:00:19.txt

echo "Start copying from SUPERFLOAT"
for filename in $( cat ${FILE_INPUT} | cut -d "," -f1 ) ; do
 do
   V1=${filename/coriolis\//}
   V2=${V1/profiles\//}
   cp -v ${ONLINE_REPO}/SUPERFLOAT/${V2} ${ONLINE_REPO}/PPCON/${V2}
 done

echo "... done"


my_prex_or_die "python -u clustering/clustering.py -i $ONLINE_REPO -u $FILE_INPUT -o $ONLINE_REPO_CLUSTERING"

my_prex_or_die "python -u make_generated_ds/generate_netcdf_netcdf4.py -t $ONLINE_REPO_CLUSTERING -m $TRAIN_DIR -p $ONLINE_REPO/PPCON"


my_prex_or_die "python dump_index.py -i $ONLINE_REPO${PPCON} -o $ONLINE_REPO${PPCON}/Float_Index.txt -t ppcon_float"

