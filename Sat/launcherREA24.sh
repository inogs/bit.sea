CHECKEDDIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/DAILY/CHECKEDfrom1999/
WEEKLY_DIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweighted
WEEKLYDATES=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYDATES
MAPCOUNT_DIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYmapcount
STD=3.5

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES
mkdir -p $MAPCOUNT_DIR

echo mpirun -np 32 python aveSatWeighted.py -i $CHECKEDDIR -o $WEEKLY_DIR -s $STD -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday -p $MAPCOUNT_DIR
