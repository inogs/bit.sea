MASKFILE=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc
INCLIMA=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_Climatology_24.nc
OUTDIR=/gpfs/scratch/userexternal/ateruzzi/KdProduct/CLIMA_FILLED/

echo python fill_climatology.py -m $MASKFILE -i $INCLIMA -o $OUTDIR
