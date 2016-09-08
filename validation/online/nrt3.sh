#! /bin/bash
MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc
LOC=/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/wrkdir/2/POSTPROC/AVE_FREQ_1/online_validation
INPUTDIR=${LOC}/PREVIOUS/TMP/
PREV_DIR=${LOC}/PREVIOUS/PROFILATORE/PROFILES
ACTU_DIR=${LOC}/ACTUAL/PROFILATORE/PROFILES

BASEDIR=NRT3
PROFILES_DIR=${BASEDIR}/PROFILES
mkdir -p $PROFILES_DIR
rm -f $PROFILES_DIR/*nc
cp $PREV_DIR/*nc $PROFILES_DIR
cp $ACTU_DIR/*nc $PROFILES_DIR #overwriting old one

python nrt3.py -o pippo.nc -b $BASEDIR -m $MASKFILE -i $INPUTDIR -d 20160628

