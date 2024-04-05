#! /bin/bash

#SBATCH --job-name=generate_superfloat
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_interactive
#SBATCH --gres=gpu:1

source ../../Sat/profile.inc
source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence3.sh


export ONLINE_REPO=/g100_scratch/userexternal/gbolzon0/carolina/ONLINE_REPO/
STARTTIME=20240209
ENDTIME=20240210
export MASKFILE=/g100_work/OGS_prodC/OPA/V10C/prod/wrkdir/analysis/2/MODEL/meshmask.nc


VALIDATION_DIR=/g100_scratch/userexternal/gbolzon0/carolina/VALIDATION

PROFILERDIR_AN=${VALIDATION_DIR}/AN_FOR_FLOATS/PROFILATORE
PROFILERDIR_FC=${VALIDATION_DIR}/FC_FOR_FLOATS/PROFILATORE
    IMG_DIR_AN=${VALIDATION_DIR}/OUTPUTS_FLOATS/ANALYSIS
    IMG_DIR_FC=${VALIDATION_DIR}/OUTPUTS_FLOATS/FORECAST
mkdir -p $IMG_DIR_AN $IMG_DIR_FC $PROFILERDIR_AN $PROFILERDIR_FC


AVEDIR_FC=/g100_scratch/userexternal/gbolzon0/carolina/FORECAST/
AVEDIR_AN=/g100_scratch/userexternal/gbolzon0/carolina/ANALYSIS/


my_prex_or_die "python -u  ../../validation/online/float_extractor_ppcon.py -st ${STARTTIME} -et ${ENDTIME} -i ${AVEDIR_AN} -b $PROFILERDIR_AN -o $IMG_DIR_AN -m $MASKFILE"
my_prex_or_die "python -u  ../../validation/online/float_extractor_ppcon.py -st ${STARTTIME} -et ${ENDTIME} -i ${AVEDIR_FC} -b $PROFILERDIR_FC -o $IMG_DIR_FC -m $MASKFILE"


# plot figures
OUTDIR=${VALIDATION_DIR}/FLOAT_MODEL_4_PROFILES
mkdir -p $OUTDIR

my_prex_or_die "python ../../validation/online/profileplotter_4.py -p ${IMG_DIR_FC}/ -a ${IMG_DIR_AN}/ -o $OUTDIR -f ${VALIDATION_DIR}/BioFloats_Descriptor.xml"

