CHECKEDDIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/DAILY/CHECKEDfrom1999/
ORIGDIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/DAILY/ORIGfrom1999/
WEEKLYbase=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLY_ORIGweight/
WEEKLY_DIR=$WEEKLYbase/WEEKLY/
WEEKLYDATES=$WEEKLYbase/WEEKLYDATES
MAPCOUNT_DIR=$WEEKLYbase/WEEKLYmapcount
STD=3.5

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES
mkdir -p $MAPCOUNT_DIR

echo mpirun -np 32 python aveSatWeighted.py -i $ORIGDIR -o $WEEKLY_DIR -s $STD -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday -p $MAPCOUNT_DIR


WEEKLYDIR_24=$WEEKLYbase/WEEKLY_24/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc
mkdir -p $WEEKLYDIR_24

echo mpirun -np 32 python interpolator.py -i $WEEKLY_DIR -o $WEEKLYDIR_24 -m Mesh24 -M $MASKFILE --inmesh SAT1km_mesh


##############################
echo NON WEIGHTED

WEEKLYbase=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYnoweight/
WEEKLY_DIR=$WEEKLYbase/WEEKLY/
WEEKLYDATES=$WEEKLYbase/WEEKLYDATES

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES

echo mpirun -np 32 python aveSat.py -i $CHECKEDDIR -o $WEEKLY_DIR -d $WEEKLYDATES -m SAT1km_mesh -t weekly_monday


WEEKLYDIR_24=$WEEKLYbase/WEEKLY_24/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc
mkdir -p $WEEKLYDIR_24

echo mpirun -np 32 python interpolator.py -i $WEEKLY_DIR -o $WEEKLYDIR_24 -m Mesh24 -M $MASKFILE --inmesh SAT1km_mesh

##############################
echo 10 DAYS AVERAGE

WEEKLYbase=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/TENDAYS/
WEEKLY_DIR=$WEEKLYbase/TENDAY/
WEEKLYDATES=$WEEKLYbase/TENDAYDATES

mkdir -p $WEEKLY_DIR
mkdir -p $WEEKLYDATES

echo mpirun -np 32 python aveSat.py -i $CHECKEDDIR -o $WEEKLY_DIR -d $WEEKLYDATES -m SAT1km_mesh -t tendays


WEEKLYDIR_24=$WEEKLYbase/TENDAY_24/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc
mkdir -p $WEEKLYDIR_24

echo mpirun -np 32 python interpolator.py -i $WEEKLY_DIR -o $WEEKLYDIR_24 -m Mesh24 -M $MASKFILE --inmesh SAT1km_mesh

##############################


echo ##############
echo verify sat variance

INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/WEEKLY/
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/WEEKLY_CHECKEDw/
INSAT=$CHECKEDDIR
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/DAILY_CHECKED/
mesh=SAT1km_mesh
INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/TESTS_forlowvarsat/SAT_10days_mario/
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/TENDAYSmario/
mesh=V4mesh
INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/WEEKLY_24/
OUTDIRP=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/CHECKpoints/WEEKLY_C24/
mesh=Mesh24
mkdir -p $OUTDIRP


echo python sat_extractpoints.py -i $INSAT -o $OUTDIRP -m $mesh -n 6 25.4 -t 41.5 33.5

