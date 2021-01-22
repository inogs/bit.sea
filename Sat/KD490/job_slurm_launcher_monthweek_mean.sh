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


CHECKED=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/KD490_checked/
OUTMONTHLY=$BASEDIR/MONTHLYAVE/
TTOUTMONTHLY=$BASEDIR/TTMONTHLYAVE/
OUTWEEKLY=$BASEDIR/WEEKLYAVE/
SATMESH=SAT1km_mesh

mkdir -p $OUTMONTHLY
mkdir -p $TTOUTMONTHLY
mkdir -p $OUTWEEKLY


my_prex_or_die "mpirun python aveSat_kd.py -i $CHECKED -o $OUTMONTHLY -m $SATMESH -t monthly"

my_prex_or_die "mpirun python aveSat_kd.py -i $CHECKED -o $OUTWEEKLY -m $SATMESH -t weekly_monday"
