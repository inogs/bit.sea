
module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea

WEEKLY_DIR=/gpfs/scratch/userexternal/gcoidess/SAT/KdProduct/WEEKLYAVE_24/
# WEEKLY_DIR contains weekly kd490 obtained from KD490 CMEMS product interpolated at model resolution
MONTHLY_DIR=/gpfs/scratch/userexternal/gcoidess/SAT/KdProduct/MONTHLYAVE_24/
# MONTHLY_DIR contains monthly kd490 obtained from KD490 CMEMS product interpolated at model resolution
CLIMA_FILE=/gpfs/scratch/userexternal/gcoidess/SAT/KdProduct/NRT_V7_newmesh/CLIMA_FILLED/
# CLIMA_FILE daily climatological file obtained by KD490 CMEMS product interpolated at model resolution and gap-filled with nearest

THEMASK=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/gdept_3d/ogstm/meshmask.nc
# THEMASK is the mask at the model resolution

######
# DA QUI IN POI COMPLETARE CON CARTELLE DI OUTPUT E 
# VALORE PERCENTUALE
#######

OUTDIR=/gpfs/scratch/userexternal/gcoidess/SAT/KdProduct/NRT_V7_newmesh/
mkdir -p $OUTDIR

OUTIND=$OUTDIR/IND24
# OUTIND will contains npy files with indexes of surface mask (used in the following script)
mkdir -p $OUTIND

WEEKLYC=$OUTDIR/DAYS7_24compl
# WEEKLYC will contains gap-filled weekly maps of KD490 at model resolution (used in the following script)
mkdir -p $WEEKLYC
 
python complete_maps.py -w $WEEKLY_DIR -d $MONTHLY_DIR -c $CLIMA_FILE -o $WEEKLYC -m $THEMASK -n $OUTIND



PERC=60 # Increasing factor = 1.18
#PERC=18 # Increasing factor = 1.18
# percentage of which Kd is increased (to keep DCM depths close to known values)
# if PERC=0 Kd is not increased

OUTDIR=$OUTDIR/OUTEOFS_p${PERC}
# OUTDIR will contains gap-filled weekly maps of KD490 obtained with EOFs signal reconstruction at model resolution
mkdir -p $OUTDIR


python eofs_map_p.py -i $WEEKLYC -o $OUTDIR -p $PERC -n $OUTIND -m $THEMASK

