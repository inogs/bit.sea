#!/bin/bash

#SBATCH --job-name=KID490_check
#SBATCH -N5
#SBATCH --ntasks-per-node=10
#SBATCH --time=1:30:00
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


date
. ./profile.inc

BASEDIR=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/KdProduct/


WEEKLY_1KMDIR=$BASEDIR/WEEKLYAVE/
MONTHLY_1KMDIR=$BASEDIR/MONTHLYAVE/
MASK24=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/gdept_3d/ogstm/meshmask.nc
WEEKLY_24_DIR=$BASEDIR/WEEKLYAVE_24/
MONTHLY_24_DIR=$BASEDIR/MONTHLYAVE_24/

mkdir -p $WEEKLY_24_DIR
mkdir -p $MONTHLY_24_DIR

my_prex_or_die "mpirun python interpolator.py -i $WEEKLY_1KMDIR -o $WEEKLY_24_DIR --inmesh SAT1km_mesh -m Mesh24 -M $MASK24" 

my_prex_or_die "mpirun python interpolator.py -i $MONTHLY_1KMDIR -o $MONTHLY_24_DIR --inmesh SAT1km_mesh -m Mesh24 -M $MASK24" 



