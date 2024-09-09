#!/bin/bash

#SBATCH --job-name=KID490_check
#SBATCH -N5
#SBATCH --ntasks-per-node=3
#SBATCH --time=03:00:00
#SBATCH --mem=115gb
#SBATCH --account=OGS_prod_0
#SBATCH -p gll_meteo_prod
#SBATCH --dependency=afterany:8354394
#SBATCH --qos=gll_qos_meteoogs

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


ORIGDIR=/gpfs/scratch/userexternal/gcoidess/SAT/DAILY_KD490/
CHECKDIR=/gpfs/scratch/userexternal/gcoidess/SAT/CHECKED/
CLIMFILE=/gpfs/scratch/userexternal/gcoidess/SAT/KD490_Climatology_1km.nc
SATMESH=SAT1km_mesh

my_prex_or_die "mpirun python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIMFILE -m $SATMESH"


date
