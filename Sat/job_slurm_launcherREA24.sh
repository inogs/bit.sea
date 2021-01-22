#!/bin/bash

#SBATCH --job-name=chl_int
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


BASEDIR=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_product/

CHECKEDDIR=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_checked/CHECKED/
#ORIGDIR=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_orig/
WEEKLYbase=$BASEDIR/WEEKLY_ORIGweight/
WEEKLY_DIR=$WEEKLYbase/WEEKLY/
WEEKLYDATES=$WEEKLYbase/WEEKLYDATES
MAPCOUNT_DIR=$WEEKLYbase/WEEKLYmapcount
STD=3.5

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES
mkdir -p $MAPCOUNT_DIR

my_prex_or_die "mpirun python aveSatWeighted.py -i $CHECKEDDIR -o $WEEKLY_DIR -s $STD -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday -p $MAPCOUNT_DIR"


WEEKLYDIR_24=$WEEKLYbase/WEEKLY_24/
MASKFILE=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/gdept_3d/ogstm/meshmask.nc
mkdir -p $WEEKLYDIR_24

my_prex_or_die "mpirun python interpolator.py -i $WEEKLY_DIR -o $WEEKLYDIR_24 -m Mesh24 -M $MASKFILE --inmesh SAT1km_mesh"

exit 0
#############################


echo ##############
echo verify sat variance

INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/WEEKLY/
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/WEEKLY_CHECKEDw/
INSAT=$CHECKEDDIR
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/DAILY_CHECKED/
mesh=SAT1km_mesh
INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/TESTS_forlowvarsat/SAT_10days_mario/
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/TENDAYSmario/
mesh=V4mesh
INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/WEEKLY_24/
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/WEEKLY_C24/
mesh=Mesh24
mkdir -p $OUTDIRP


echo python sat_extractpoints.py -i $INSAT -o $OUTDIRP -m $mesh -n 6 25.4 -t 41.5 33.5

