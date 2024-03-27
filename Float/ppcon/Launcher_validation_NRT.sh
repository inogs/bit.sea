#! /bin/bash

#SBATCH --job-name=generate_superfloat
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_interactive
#SBATCH --gres=gpu:1

source /g100_scratch/userexternal/camadio0/PPCON/bit.sea/Sat/profile.inc
my_prex_or_die "command"

#/g100_work/OGS_prodC/OPA/V10C/prod/inpdir/ONLINE/SUPERFLOAT/
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18/bin/activate
ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/
STARTTIME="20240209"
ENDTIME="20240210"
export ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/PPCON/

############################################################
################ 1. parsing per validation #################
############################################################

echo "The dataset is ready to be used..."
date
echo "..start validation .. "

export PYTHONPATH=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/:$PYTHONPATH
export ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/PPCON/
export MASKFILE=/g100_work/OGS_prodC/OPA/V10C/prod/wrkdir/analysis/2/MODEL/meshmask.nc
BASEPPCON=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/validation/online/

mkdir -p ${BASEPPCON}FC_FOR_FLOATS/
mkdir -p ${BASEPPCON}OUTPUTS_FLOATS/
############################################################
##### 2. directory I/O forecast and extraction [mathcup] ###
############################################################
BIO_DIRP=/g100_scratch/userexternal/gbolzon0/carolina/FORECAST/
PROFILERDIRP=${BASEPPCON}/FC_FOR_FLOATS/PROFILATORE
IMG_DIR_FOR=${BASEPPCON}/OUTPUTS_FLOATS/FORECAST
mkdir -p $IMG_DIR_FOR $PROFILERDIRP
MASKFILE=/g100_work/OGS_prodC/OPA/V10C/prod/wrkdir/analysis/2/MODEL/meshmask.nc
#my_prex_or_die "python -u ${BASEPPCON}/float_extractor_ppcon.py -st ${STARTTIME} -et ${ENDTIME} -i ${BIO_DIRP} -b $PROFILERDIRP -o $IMG_DIR_FOR -m $MASKFILE"


# directory I/O analysis and extraction [mathcup]
BIO_DIRP=/g100_scratch/userexternal/gbolzon0/carolina/ANALYSIS/
PROFILERDIRP=${BASEPPCON}/ANALISYS/PROFILATORE
IMG_DIR_AN=${BASEPPCON}/OUTPUTS_FLOATS/ANALISYS/
mkdir -p ${BASEPPCON}/ANALISYS/
mkdir -p $IMG_DIR_AN $PROFILERDIRP
#my_prex_or_die "python -u ${BASEPPCON}/float_extractor_ppcon.py -st ${STARTTIME} -et ${ENDTIME} -i ${BIO_DIRP} -b $PROFILERDIRP -o $IMG_DIR_AN -m $MASKFILE"

# plot figures
OUTIR=${BASEPPCON}/MEDEAF_FLOAT_MODEL_PROFILE
mkdir -p $OUTIR

my_prex_or_die "python ${BASEPPCON}/profileplotter_4.py -p ${IMG_DIR_FOR}/ -a ${IMG_DIR_AN}/ -o $OUTIR -f ${BASEPPCON}/BioFloats_Descriptor.xml"

date
