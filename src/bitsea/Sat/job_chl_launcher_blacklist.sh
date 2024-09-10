#!/bin/bash

#SBATCH --job-name=chl_check
#SBATCH -N5
#SBATCH --ntasks-per-node=3
#SBATCH --time=1:30:00
#SBATCH --mem=115gb
#SBATCH --account=OGS20_PRACE_P
#SBATCH --partition=gll_usr_prod

#cd $SLURM_SUBMIT_DIR

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea


. ./profile.inc

date


#annaDIR=/pico/scratch/userexternal/ateruzzi/Sat_Check_Clima/

  DIR_ORIG1km=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_orig/
DIR_CHECK_1km=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/chl_checked/
#DIR_CLIMA1km=
   FILEClima=/gpfs/scratch/userexternal/gcoidess/SAT/bit.sea/Sat/Climatology_CHL.nc
 DIR_SUBMASK=/gpfs/scratch/userexternal/gcoidess/TMP_DOWNLOAD/SUBMASKsat/

#DIR_CHECK_1km=/marconi_scratch/userexternal/gbolzon0/SAT_CHECK

## Creation of subask indexes for sat mesh
# To be executed once and for all (once for each climatology)
my_prex_or_die "python sat_indsub.py -o $DIR_SUBMASK -c $FILEClima -m SAT1km_mesh -l /gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/gdept_3d/ogstm/meshmask.nc -t All"
my_prex_or_die "python sat_indsub.py -o $DIR_SUBMASK -c $FILEClima -m SAT1km_mesh -l /gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/gdept_3d/ogstm/meshmask.nc -t Open"
# --------------------------------
## Creation of climatology statistics for the Mediterranean subbasins
# To be executed once and for all (once for each climatology)
# OUT_STATSCLIM=$annaDIR/STATS_CLIM/
# mkdir $OUT_STATSCLIM
# echo python stats_sub_clima.py -c $FILEClima -m SAT1km_mesh  -s $DIR_SUBMASK -w $OUT_STATSCLIM

# --------------------------------
## Check of sat files

OUT_STATS=$DIR_CHECK_1km/STATISTICS
mkdir -p $DIR_CHECK_1km/CHECKED
mkdir -p $DIR_CHECK_1km/REJECTED
mkdir -p $OUT_STATS

my_prex_or_die "mpirun python sat_check.py -i $DIR_ORIG1km -o $DIR_CHECK_1km -c $FILEClima -m SAT1km_mesh  -s $DIR_SUBMASK -w $OUT_STATS"


exit 0
# --------------------------------
## Figures

# OUT_FIG=FIGURES/$TYPEIN/
# mkdir -p $OUT_FIG
# my_prex_or_die "mpirun python plot_climastats_time.py -o $OUT_FIG -i $OUT_STATS -m $DIR_SUBMASK"


# --------------------------------
## Ave sat

WEEKLYMULTI=$annaDIR/WEEKLY
WEEKLYDATES=$annaDIR/WEEKLY_2_AVEDATES
mkdir -p $WEEKLYMULTI
mkdir -p $WEEKLYDATES

my_prex_or_die "mpirun python aveSat.py -i $DIR_CHECK_1km/CHECKED -o $WEEKLYMULTI -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday"


# --------------------------------
## Compose blacklisting file

FLAGSATDIR=$annaDIR/FakeArchive
OUTNC=$annaDIR/NC_OUT
MASKMOD=/pico/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
mkdir -p $OUTNC

my_prex_or_die "mpirun python blacklisting_ncfile.py -r $DIR_CHECK_1km/REJECTED -d $FLAGSATDIR -s $WEEKLYDATES -o $OUTNC -t SAT1km_mesh -m $MASKMOD"

