
WEEKLY_DIR=/gpfs/scratch/userexternal/ateruzzi/KdProduct/WEEKLYAVE_24/
# WEEKLY_DIR contains weekly kd490 obtained from KD490 CMEMS product interpolated at model resolution
MONTHLY_DIR=/gpfs/scratch/userexternal/ateruzzi/KdProduct/MONTHLYAVE_24/
# MONTHLY_DIR contains monthly kd490 obtained from KD490 CMEMS product interpolated at model resolution
CLIMA_FILE=/gpfs/scratch/userexternal/ateruzzi/KdProduct/NRT_V7_newmesh/CLIMA_FILLED/KD490_Climatology_24_filled.nc
# CLIMA_FILE daily climatological file obtained by KD490 CMEMS product interpolated at model resolution and gap-filled with nearest

THEMASK=/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask.nc
# THEMASK is the mask at the model resolution

######
# DA QUI IN POI COMPLETARE CON CARTELLE DI OUTPUT E 
# VALORE PERCENTUALE
#######

OUTDIR=/gpfs/scratch/userexternal/ateruzzi/KdProduct/NRT_V7_newmesh/
OUTIND=$OUTDIR/IND24
# OUTIND will contains npy files with indexes of surface mask (used in the following script)
mkdir -p $OUTIND

WEEKLYC=$OUTDIR/DAYS7_24compl
# WEEKLYC will contains gap-filled weekly maps of KD490 at model resolution (used in the following script)
mkdir -p $WEEKLYC
 
echo python complete_maps.py -w $WEEKLY_DIR -d $MONTHLY_DIR -c $CLIMA_FILE -o $WEEKLYC -m $THEMASK -n $OUTIND


PERC=60 # Increasing factor = 1.18
#PERC=18 # Increasing factor = 1.18
# percentage of which Kd is increased (to keep DCM depths close to known values)
# if PERC=0 Kd is not increased

OUTDIR=$OUTDIR/OUTEOFS_p${PERC}
# OUTDIR will contains gap-filled weekly maps of KD490 obtained with EOFs signal reconstruction at model resolution
mkdir -p $OUTDIR


echo python eofs_map_p.py -i $WEEKLYC -o $OUTDIR -p $PERC -n $OUTIND -m $THEMASK

