MASKFILE=/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask.nc
INCLIMA=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/KD490_Climatology_24.nc
OUTDIR=/gpfs/scratch/userexternal/ateruzzi/KdProduct/NRT_V7_newmesh/CLIMA_FILLED/
mkdir -p $OUTDIR

echo python fill_climatology.py -m $MASKFILE -i $INCLIMA -o $OUTDIR
