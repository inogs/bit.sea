#!/bin/bash

#SBATCH --job-name=POST
#SBATCH -N5
#SBATCH --ntasks-per-node=10
#SBATCH --time=01:30:00
#SBATCH --mem=115gb
#SBATCH --account=OGS20_PRACE_P
#SBATCH --partition=gll_usr_prod

cd $SLURM_SUBMIT_DIR

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea


. ./profile.inc


BASEDIR=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_product/

CHECKED_DIR=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_checked/CHECKED/
MONTHLY_DIR=$BASEDIR/MONTHLY_2019
DIR_DATE=$BASEDIR/MONTHLY_dates_2019

mkdir -p $MONTHLY_DIR
mkdir -p $DIR_DATE

my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR -o $MONTHLY_DIR -m SAT1km_mesh -t monthly -d $DIR_DATE"

WEEKLY_DIR=$BASEDIR/WEEKLY_2019
DIR_DATE=$BASEDIR/WEEKLY_dates_2019

mkdir -p $WEEKLY_DIR
mkdir -p $DIR_DATE

my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR -o $MONTHLY_DIR -m SAT1km_mesh -t weekly_monday -d $DIR_DATE"




