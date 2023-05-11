#! /bin/bash

#SBATCH --job-name=CHL_all
#SBATCH -N1
#SBATCH --ntasks-per-node=12
#SBATCH --time=01:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_devC
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo


cd $SLURM_SUBMIT_DIR

. ./profile.inc

source /g100_work/OGS21_PRACE_P/COPERNICUS/sequence3.sh

date

BASEDIR=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/INPUT_FILES/CHL/

#ORIGDIR=$BASEDIR/ORIG/
ORIGDIR=$BASEDIR/ORIG_del_mese
CHECKDIR=$BASEDIR/DAILY/CHECKED/
REJECTED_DIR=$BASEDIR/DAILY/REJECTED/
CHECK_BASEDIR=$BASEDIR/DAILY/
STATDIR=$BASEDIR/DAILY/STATISTICS/
WEEKLY_DIR1km=$BASEDIR/WEEKLY_1_1km/
WEEKLYDATES=$BASEDIR/WEEKLY_1_AVEDATES/
WEEKLY_DIR_24=$BASEDIR/WEEKLY_1_24/
CLIM_FILE=$BASEDIR/../STATIC/SatClimatology.nc
DIR_SUBMASK=$BASEDIR/../STATIC/SUBMASKsat/
MASKFILE=$BASEDIR/../STATIC/meshmask.nc
mkdir -p $CHECKDIR
mkdir -p $CHECK_BASEDIR
mkdir -p $REJECTED_DIR
mkdir -p $WEEKLY_DIR1km
mkdir -p $WEEKLYDATES
mkdir -p $WEEKLY_DIR_24
mkdir -p $STATDIR

my_prex_or_die "mpirun python sat_check.py -i $ORIGDIR -o $CHECK_BASEDIR -c $CLIM_FILE -m SAT1km_mesh -s $DIR_SUBMASK -w $STATDIR "

my_prex_or_die "mpirun python aveSat.py -i $CHECKDIR -o $WEEKLY_DIR1km -m SAT1km_mesh -t weekly_monday -d $WEEKLYDATES "

my_prex_or_die "mpirun python interpolator.py -i $WEEKLY_DIR1km -o $WEEKLY_DIR_24 -m $MASKFILE --inmesh SAT1km_mesh "

