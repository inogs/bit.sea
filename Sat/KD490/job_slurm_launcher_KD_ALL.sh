#! /bin/bash

#SBATCH --job-name=KD_MEAN
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
. ./profile.inc

BASEDIR=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/INPUT_FILES/KD490/

mkdir -p $BASEDIR


#ORIGDIR=$BASEDIR/ORIG/
ORIGDIR=$BASEDIR/ORIG_del_mese/
CHECKDIR=$BASEDIR/CHECKED/
CLIMFILE=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/BITSEA/Sat/KD490/KD490_Climatology_1km.nc
SATMESH=SAT1km_mesh

mkdir -p $CHECKDIR

my_prex_or_die "mpirun python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIMFILE -m $SATMESH"


CHECKED=$BASEDIR/CHECKED/
OUTMONTHLY=$BASEDIR/MONTHLYAVE/
TTOUTMONTHLY=$BASEDIR/TTMONTHLYAVE/
OUTWEEKLY=$BASEDIR/WEEKLYAVE/
SATMESH=SAT1km_mesh

mkdir -p $OUTMONTHLY
mkdir -p $TTOUTMONTHLY
mkdir -p $OUTWEEKLY


my_prex_or_die "mpirun python aveSat_kd.py -i $CHECKED -o $OUTMONTHLY -m $SATMESH -t monthly"

my_prex_or_die "mpirun python aveSat_kd.py -i $CHECKED -o $OUTWEEKLY -m $SATMESH -t weekly_monday"

date

WEEKLY_1KMDIR=$BASEDIR/WEEKLYAVE/
MONTHLY_1KMDIR=$BASEDIR/MONTHLYAVE/
MASK24=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/INPUT_FILES/STATIC/meshmask.nc
WEEKLY_24_DIR=$BASEDIR/WEEKLYAVE_24/
MONTHLY_24_DIR=$BASEDIR/MONTHLYAVE_24/

mkdir -p $WEEKLY_24_DIR
mkdir -p $MONTHLY_24_DIR

my_prex_or_die "mpirun python interpolator.py -i $WEEKLY_1KMDIR -o $WEEKLY_24_DIR --inmesh SAT1km_mesh -m Mesh24 -M $MASK24"

my_prex_or_die "mpirun python interpolator.py -i $MONTHLY_1KMDIR -o $MONTHLY_24_DIR --inmesh SAT1km_mesh -m Mesh24 -M $MASK24"

date


WEEKLY_DIR=$BASEDIR/WEEKLYAVE_24/
# WEEKLY_DIR contains weekly kd490 obtained from KD490 CMEMS product interpolated at model resolution
MONTHLY_DIR=$BASEDIR/MONTHLYAVE_24/
# MONTHLY_DIR contains monthly kd490 obtained from KD490 CMEMS product interpolated at model resolution
CLIMA_FILE=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/INPUT_FILES/KD490/NRT_V7_newmesh/CLIMA_FILLED/

# CLIMA_FILE daily climatological file obtained by KD490 CMEMS product interpolated at model resolution and gap-filled with nearest
#mkdir -p $CLIMA_FILE

THEMASK=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/INPUT_FILES/STATIC/meshmask.nc
# THEMASK is the mask at the model resolution

######
# DA QUI IN POI COMPLETARE CON CARTELLE DI OUTPUT E
# VALORE PERCENTUALE
#######

OUTDIR=$BASEDIR/NRT_V7_newmesh/
mkdir -p $OUTDIR

OUTIND=$OUTDIR/IND24
# OUTIND will contains npy files with indexes of surface mask (used in the following script)
mkdir -p $OUTIND

WEEKLYC=$OUTDIR/DAYS7_24compl
# WEEKLYC will contains gap-filled weekly maps of KD490 at model resolution (used in the following script)
mkdir -p $WEEKLYC

my_prex_or_die "python complete_maps.py -w $WEEKLY_DIR -d $MONTHLY_DIR -c $CLIMA_FILE -o $WEEKLYC -m $THEMASK -n $OUTIND"

PERC=60 # Increasing factor = 1.18
#PERC=18 # Increasing factor = 1.18
# percentage of which Kd is increased (to keep DCM depths close to known values)
# if PERC=0 Kd is not increased

OUTDIR=$OUTDIR/OUTEOFS_p${PERC}
# OUTDIR will contains gap-filled weekly maps of KD490 obtained with EOFs signal reconstruction at model resolution
mkdir -p $OUTDIR


my_prex_or_die "python eofs_map_p.py -i $WEEKLYC -o $OUTDIR -p $PERC -n $OUTIND -m $THEMASK"

date
