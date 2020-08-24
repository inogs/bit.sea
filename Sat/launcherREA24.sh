CHECKEDDIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/DAILY/CHECKEDfrom1999/
WEEKLYbase=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/
WEEKLY_DIR=$WEEKLYbase/WEEKLY/
WEEKLYDATES=$WEEKLYbase/WEEKLYDATES
MAPCOUNT_DIR=$WEEKLYbase/WEEKLYmapcount
STD=3.5

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES
mkdir -p $MAPCOUNT_DIR

echo mpirun -np 32 python aveSatWeighted.py -i $CHECKEDDIR -o $WEEKLY_DIR -s $STD -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday -p $MAPCOUNT_DIR


WEEKLYDIR_24=$WEEKLYbase/WEEKLY_24/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc
mkdir -p $WEEKLYDIR_24

echo mpirun -np 32 python interpolator.py -i $WEEKLY_DIR -o $WEEKLYDIR_24 -m Mesh24 -M $MASKFILE --inmesh SAT1km_mesh



echo NON WEIGHTED

WEEKLYbase=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYnoweight/
WEEKLY_DIR=$WEEKLYbase/WEEKLY/
WEEKLYDATES=$WEEKLYbase/WEEKLYDATES
MAPCOUNT_DIR=$WEEKLYbase/WEEKLYmapcount
STD=3.5

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES
mkdir -p $MAPCOUNT_DIR

echo mpirun -np 32 python aveSat.py -i $CHECKEDDIR -o $WEEKLY_DIR -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday


WEEKLYDIR_24=$WEEKLYbase/WEEKLY_24/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc
mkdir -p $WEEKLYDIR_24

echo mpirun -np 32 python interpolator.py -i $WEEKLY_DIR -o $WEEKLYDIR_24 -m Mesh24 -M $MASKFILE --inmesh SAT1km_mesh
