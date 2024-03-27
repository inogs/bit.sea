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
BASEPPCON=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/Float/ppcon/
ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/
STARTTIME="20240325"
ENDTIME="20240326"
ONLINE_REPO_CLUSTERING=${ONLINE_REPO}clustering/

############################################################
############# 1 Creo il file tensor csv ####################
#########salva in clustering/clustering..csv ###############
############################################################

DIFFDIR=/g100_work/OGS_prodC/OPA/V10C-prod/log/daily/
cp -r ${DIFFDIR}/DIFF_floats*$STARTTIME*.txt ${ONLINE_REPO}DIFF_floats_${STARTTIME}.txt
FILE_INPUT=${ONLINE_REPO}DIFF_floats_${STARTTIME}.txt    # DIFF FILE  ––> Update dataset
my_prex_or_die "python -u clustering/clustering.py -i $ONLINE_REPO -u $FILE_INPUT -o $ONLINE_REPO_CLUSTERING"


############################################################
############# 2. creo il dataset chiamato : PPCON ##########
############################################################
PPCON=PPCON/
PPCON_MODELDIR=$BASEPPCON/results/    #inputs
mkdir -p $ONLINE_REPO${PPCON}  # creo PPCON 
# --> GB:  cp -r $ONLINE_REPO${SUPERFLOAT}/* $ONLINE_REPO${PPCON} 
cp -r ${ONLINE_REPO}/SUPERFLOAT/* $ONLINE_REPO${PPCON}

my_prex_or_die "python -u make_generated_ds/generate_netcdf_netcdf4.py -t $ONLINE_REPO_CLUSTERING -m $PPCON_MODELDIR -p $ONLINE_REPO${PPCON}"
#                                                      |legge cluster.csv | legge /results | PPCON dataset


############################################################
########## 3. creo il float indexl dataset : PPCON #########
############################################################
my_prex_or_die "python dump_index.py -i $ONLINE_REPO${PPCON} -o $ONLINE_REPO${PPCON}/Float_Index.txt -t ppcon_float"


exit 0

#______ dataset is done and vars corrected as fo r coriolis ____#

############################################################
##### 4. dovrei correggere i vals come in SUPERFLOAT #######
############################################################
#

export ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/PPCON/

############################################################
################ 5. parsing per validation #################
############################################################

echo "The dataset is ready to be used..."
date
echo "..start validation .. "

export PYTHONPATH=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/:$PYTHONPATH
export ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/PPCON/
export MASKFILE=/g100_work/OGS_prodC/OPA/V10C/prod/wrkdir/analysis/2/MODEL/meshmask.nc
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18/bin/activate

BASEPPCON=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/validation/online/

mkdir -p ${BASEPPCON}FC_FOR_FLOATS/
mkdir -p ${BASEPPCON}OUTPUTS_FLOATS/

# directory I/O forecast and extraction [mathcup]
BIO_DIRP=/g100_scratch/userexternal/gbolzon0/carolina/FORECAST/
PROFILERDIRP=${BASEPPCON}/FC_FOR_FLOATS/PROFILATORE
IMG_DIR_FOR=${BASEPPCON}/OUTPUTS_FLOATS/FORECAST
mkdir -p $IMG_DIR_FOR $PROFILERDIRP
MASKFILE=/g100_work/OGS_prodC/OPA/V10C/prod/wrkdir/analysis/2/MODEL/meshmask.nc
#python -u ${BASEPPCON}/float_extractor_ppcon.py -st ${STARTTIME} -et ${ENDTIME} -i ${BIO_DIRP} -b $PROFILERDIRP -o $IMG_DIR_FOR -m $MASKFILE

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

echo my_prex_or_die "python ${BASEPPCON}/profileplotter_4.py -p ${IMG_DIR_FOR}/ -a ${IMG_DIR_AN}/ -o $OUTIR -f ${BASEPPCON}/BioFloats_Descriptor.xml"

date
