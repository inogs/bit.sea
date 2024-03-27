#! /bin/bash

#SBATCH --job-name=generate_superfloat
#SBATCH --ntasks=1
#SBATCH --time=01:05:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_interactive
#SBATCH --gres=gpu:1

source /g100_scratch/userexternal/camadio0/PPCON/bit.sea/Sat/profile.inc
my_prex_or_die "command"

#/g100_work/OGS_prodC/OPA/V10C/prod/inpdir/ONLINE/SUPERFLOAT/
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18/bin/activate
BASEPPCON=/g100_scratch/userexternal/camadio0/PPCON/bit.sea/Float/ppcon/
ONLINE_REPO=/g100_scratch/userexternal/camadio0/ONLINE_REPO_202403/
ONLINE_REPO_CLUSTERING=${ONLINE_REPO}clustering/



############################################################
############# 1 Creo il file tensor csv ####################
###############salva in dir cluster ########################
############################################################

FILE_INPUT=${ONLINE_REPO}SUPERFLOAT/Float_Index.txt   
my_prex_or_die "python -u clustering/clustering.py -i $ONLINE_REPO -u $FILE_INPUT -o $ONLINE_REPO_CLUSTERING"





############################################################
############# 2. creo il dataset chiamato : PPCON ##########
############################################################
PPCON=PPCON/
PPCON_MODELDIR=$BASEPPCON/results/    #inputs
mkdir -p $ONLINE_REPO${NEW_DIR}
# --> GB:  cp -r $ONLINE_REPO${SUPERFLOAT}/* $ONLINE_REPO${PPCON}
cp -r ${ONLINE_REPO}/SUPERFLOAT/* $ONLINE_REPO${PPCON}

my_prex_or_die "python -u make_generated_ds/generate_netcdf_netcdf4.py -t $ONLINE_REPO_CLUSTERING -m $PPCON_MODELDIR -p $ONLINE_REPO${PPCON}"
#                                                      |legge cluster.csv | legge /results | PPCON dataset


############################################################
########## 3. creo il float indexl dataset : PPCON #########
############################################################
my_prex_or_die "python dump_index.py -i $ONLINE_REPO${PPCON} -o $ONLINE_REPO${PPCON}/Float_Index.txt -t ppcon_float"


exit 0
