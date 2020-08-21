ORIGDIR= # repository of downloaded kd490
CHECKDIR= # output dir
CLIMFILE=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_Climatology_1km.nc
SATMESH= # The same used as "-m argument" for sat_check on chlorophyll

echo python sat_check_kd.py -i $ORIGDIR -o $CHECKDIR -c $CLIMFILE -m $SATMESH
