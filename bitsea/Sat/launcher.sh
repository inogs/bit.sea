#! /bin/bash
#SBATCH --job-name=SAT
#SBATCH -N1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_devC
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo

cd $SLURM_SUBMIT_DIR

. ./profile.inc

source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence3.sh

unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false
MASKFILE=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc

export ONLINE_REPO=/g100_scratch/userexternal/gbolzon0/V10C


CHECKED_DIR=${ONLINE_REPO}/SAT/CHL/DT/DAILY/CHECKED
CHECKED_DIR24=${ONLINE_REPO}/SAT/CHL/DT/DAILY/CHECKED_24
WEEKLY_DIR1km=${ONLINE_REPO}/SAT/CHL/DT/WEEKLY_1_1km
WEEKLY_DIR_24=${ONLINE_REPO}/SAT/CHL/DT/WEEKLY_1_24/
WEEKLY_DIR424=${ONLINE_REPO}/SAT/CHL/DT/WEEKLY_4_24

mkdir -p $WEEKLY_DIR1km $WEEKLY_DIR_24 $WEEKLY_DIR424
my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR1km -m SAT1km_mesh -t weekly_monday -v CHL"
my_prex_or_die "mpirun python interpolator.py -i $WEEKLY_DIR1km -o $WEEKLY_DIR_24 -m $MASKFILE --inmesh SAT1km_mesh -v CHL"

for pftvar in CHL DIATO NANO PICO DINO ; do
   my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR24 -o $WEEKLY_DIR424 -m Mesh24 -t weekly_thursday -v $pftvar"
done


CHECKED_DIR=${ONLINE_REPO}/SAT/KD490/DT/DAILY/CHECKED
CHECKED_DIR24=${ONLINE_REPO}/SAT/KD490/DT/DAILY/CHECKED_24
WEEKLY_DIR1km=${ONLINE_REPO}/SAT/KD490/DT/WEEKLY_1_1km
WEEKLY_DIR_24=${ONLINE_REPO}/SAT/KD490/DT/WEEKLY_1_24/
WEEKLY_DIR424=${ONLINE_REPO}/SAT/KD490/DT/WEEKLY_4_24

mkdir -p $WEEKLY_DIR1km $WEEKLY_DIR_24 $WEEKLY_DIR424
my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR1km -m SAT1km_mesh -t weekly_monday -v KD490"
my_prex_or_die "mpirun python interpolator.py -i $WEEKLY_DIR1km -o $WEEKLY_DIR_24 -m $MASKFILE --inmesh SAT1km_mesh -v KD490"
my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR24 -o $WEEKLY_DIR424 -m Mesh24 -t weekly_thursday -v KD490"

CHECKED_DIR24=${ONLINE_REPO}/SAT/KD490/NRT/DAILY/CHECKED_24
WEEKLY_DIR424=${ONLINE_REPO}/SAT/KD490/NRT/WEEKLY_4_24
mkdir -p $WEEKLY_DIR424
my_prex_or_die "mpirun python aveSat.py -i $CHECKED_DIR24 -o $WEEKLY_DIR424 -m Mesh24 -t weekly_thursday -v KD490"




exit 0


WEEKLY_DIR_24=/g100_scratch/userexternal/gbolzon0/V10C/SAT/CHL/DT/WEEKLY_4_24
mkdir -p $WEEKLY_DIR_24
for pftvar in CHL DIATO NANO PICO DINO ; do
   my_prex_or_die "mpirun -np 24 python interpolator.py -i $WEEKLY_DIR1km -o $WEEKLY_DIR_24 -m $MASKFILE --inmesh SAT1km_mesh -v $pftvar "
done

exit 0






#####  REFLECTANCE SECTION ###########
ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/OCEANCOLOUR_MED_BGC_L3_MY_009_143/preproc_download/Download/my.cmems-du.eu/Core/OCEANCOLOUR_MED_BGC_L3_MY_009_143/cmems_obs-oc_med_bgc-reflectance_my_l3-multi-1km_P1D/2019/01
       ORIGDIR=/g100_scratch/userexternal/gbolzon0/V11C/SAT/DAILY/ORIG
   CHECKED_DIR=/g100_scratch/userexternal/gbolzon0/V11C/SAT/DAILY/CHECKED
 WEEKLY_DIR1km=/g100_scratch/userexternal/gbolzon0/V11C/SAT/WEEKLY_1km
 WEEKLY_DIR_24=/g100_scratch/userexternal/gbolzon0/V11C/SAT/WEEKLY_24


mkdir -p $CHECKED_DIR $WEEKLY_DIR1km $WEEKLY_DIR_24
MPI="mpirun -np 24"
for nm in 412 443 490 510 555 670; do
   var=RRS${nm}
   my_prex_or_die "$MPI python sat_check_QI.py -i $ORIGDIR -o $CHECKED_DIR -m SAT1km_mesh -v $var --QI 3.0"
   my_prex_or_die "$MPI python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR1km -m SAT1km_mesh -t weekly_thursday -v $var"
   my_prex_or_die "$MPI python interpolator.py -i $WEEKLY_DIR1km -o $WEEKLY_DIR_24 -m $MASKFILE --inmesh SAT1km_mesh -v $var "

done

exit 0

##############   KD490 SECTION ##########################################################
ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/SAT/KD490/DT/DAILY/ORIG/
CHECKED_DIR=/g100_scratch/userexternal/gbolzon0/V10C/SAT/KD490/DT/DAILY/CHECKED/SIGMA_2.0/

mkdir -p $CHECKED_DIR
my_prex_or_die "python sat_check_QI.py -i $ORIGDIR -o $CHECKED_DIR -m SAT1km_mesh -v KD490 --QI 2.0 --Kd_min 0.021"

WEEKLY_DIR1km=/g100_scratch/userexternal/gbolzon0/V10C/SAT/KD490/DT/SIGMA_2.0/WEEKLY_4_1km
mkdir -p $WEEKLY_DIR1km

my_prex_or_die " python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR1km -m SAT1km_mesh -t weekly_thursday -v KD490"

WEEKLY_DIR_24=/g100_scratch/userexternal/gbolzon0/V10C/SAT/KD490/DT/SIGMA_2.0/WEEKLY_4_24
mkdir -p $WEEKLY_DIR_24

my_prex_or_die " python interpolator.py -i $WEEKLY_DIR1km -o $WEEKLY_DIR_24 -m $MASKFILE --inmesh SAT1km_mesh -v KD490 "




ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/SAT/CHL/DT/DAILY/ORIG/

CHECKED_DIR=/g100_scratch/userexternal/gbolzon0/V10C/SAT/CHL/DT/DAILY/CHECKED



#my_prex_or_die "python sat_check_QI.py -i $ORIGDIR -o $CHECKED_DIR -m SAT1km_mesh -v CHL --QI 2"
#my_prex_or_die "python sat_check_QI.py -i $ORIGDIR -o $CHECKED_DIR -m SAT1km_mesh -v DIATO -f --QI 2 "

for pftvar in CHL DIATO NANO PICO DINO ; do
    break
    my_prex_or_die "python sat_check_QI.py -i $ORIGDIR -o $CHECKED_DIR -m SAT1km_mesh -v $pftvar --QI 2 "
done



