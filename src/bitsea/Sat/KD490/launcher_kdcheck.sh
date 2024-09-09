ORIGDIR=/gpfs/scratch/userexternal/gcoidess/SAT/DAILY_KD490/ # repository of downloaded kd490
CHECKDIR=/gpfs/scratch/userexternal/gcoidess/SAT/CHECKED/

mkdir -p $CHECKDIR

CLIMFILE=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_Climatology_1km.nc
SATMESH=KD490mesh # The same used as "-m argument" for sat_check on chlorophyll

python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIMFILE -m $SATMESH
