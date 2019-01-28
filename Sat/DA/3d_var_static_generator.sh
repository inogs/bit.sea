#! /bin/bash


MASKFILE=/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc
 SAT_DIR_10=/pico/scratch/userexternal/gbolzon0/EOF-python/CCI1km_Interp24_10days/
  VARSATDIR=/pico/scratch/userexternal/gbolzon0/EOF-python/VARSAT/
VARSAT_INCR=/pico/scratch/userexternal/gbolzon0/EOF-python/SUMMER_INCREASED50_100/



python varSat.py -i $SAT_DIR_10 -o $VARSATDIR -m Mesh24

python increaseSummerSatVar.py -i $VARSATDIR -o $VARSAT_INCR -m $MASKFILE


     MODELDIR=/pico/scratch/userexternal/gbolzon0/EOF-python/ORIG_INPUTS/CHL_SUP24/
SATDIR_WEEKLY=/pico/scratch/userexternal/gbolzon0/EOF-python/ORIG_INPUTS/CCI1km_Interp24_weekly/
  VAR_ERR_DIR=/pico/scratch/userexternal/gbolzon0/EOF-python/VarErr/

python varErr.py -i $MODELDIR -s $SATDIR_WEEKLY -o $VAR_ERR_DIR -m Mesh24

   VAR_MOD_DIR=/pico/scratch/userexternal/gbolzon0/EOF-python/VarMod/

python varMod.py -i $VAR_ERR_DIR -s $VARSATDIR -o $VAR_MOD_DIR -m $MASKFILE
