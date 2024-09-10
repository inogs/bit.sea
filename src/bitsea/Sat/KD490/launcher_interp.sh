WEEKLY_1KMDIR=$CINECA_SCRATCH/KdProduct/WEEKLYAVE/
MONTHLY_1KMDIR=$CINECA_SCRATCH/KdProduct/MONTHLYAVE/
MASK24=/gpfs/scratch/userexternal/gcoidess/REA_24/TEST_09/wrkdir/MODEL/meshmask.nc
WEEKLY_24_DIR=$CINECA_SCRATCH/KdProduct/WEEKLYAVE_24/
MONTHLY_24_DIR=$CINECA_SCRATCH/KdProduct/MONTHLYAVE_24/

mkdir -p $WEEKLY_24_DIR
mkdir -p $MONTHLY_24_DIR

echo mpirun -np 35 python interpolator.py -i $WEEKLY_1KMDIR -o $WEEKLY_24_DIR --inmesh SAT1km_mesh -m Mesh24 -M $MASK24 

echo mpirun -np 35 python interpolator.py -i $MONTHLY_1KMDIR -o $MONTHLY_24_DIR --inmesh SAT1km_mesh -m Mesh24 -M $MASK24 
