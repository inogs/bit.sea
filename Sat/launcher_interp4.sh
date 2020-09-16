
#DAILY
# Interploation from 1/16 to 1/4 of REP SAT 
DAILY_1km_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_4/DAILYrep1km_2013/
DAILY_4_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_4/DAILYrep4_2013/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/SAT_4/meshmask_4.nc


mkdir -p $DAILY_4_DIR
echo mpirun -np 35 python interpolator.py -i $DAILY_1km_DIR -o $DAILY_4_DIR --inmesh SAT1km_mesh  -m Mesh4 -M $MASKFILE

exit 0

#WEEKLY
# Interploation from 1/16 to 1/4 of REP SAT 
WEEKLY_16_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_16_2013/WEEKLY/
WEEKLY_4_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_4/WEEKLY/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS4/meshmask.nc


mkdir -p $WEEKLY_4_DIR
echo mpirun -np 10 python interpolator.py -i $WEEKLY_16_DIR -o $WEEKLY_4_DIR --inmesh V4mesh -m Mesh4 -M $MASKFILE




