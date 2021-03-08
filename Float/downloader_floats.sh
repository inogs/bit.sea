#!/bin/bash

DIFF_FILE=/galileo/home/userexternal/gcoidess/test_float_download/DIFF_floats.txt
LOCALDIR=/gpfs/scratch/userexternal/gbolzon0/V7C/TMP_CORIOLIS
HOST=ftp.ifremer.fr
USER=anonymous
REMOTEDIR=/ifremer/argo/dac/
export ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/V7C/ONLINE

mkdir -p $LOCALDIR
for I in `awk --field-separator=, '{print $1}'  $DIFF_FILE ` ; do
  ncftpget -u $USER $HOST $LOCALDIR $REMOTEDIR/$I
done


echo "Moving NetCDF files from TMP dir to ONLINE CORIOLIS file"
cd $LOCALDIR
for I in `ls *nc ` ; do
   WMO=${I:2:7}
   CORIOLIS_DIR=$ONLINE_REPO/CORIOLIS/${WMO}/
   mkdir -p $CORIOLIS_DIR
   mv $I $CORIOLIS_DIR
done