CHECKED=/gpfs/scratch/userexternal/gcoidess/REA_24/REA_24_SAT/KD490/DAILY/CHECKED/
OUTMONTHLY=$CINECA_SCRATCH/KdProduct/MONTHLYAVE/
TTOUTMONTHLY=$CINECA_SCRATCH/KdProduct/TTMONTHLYAVE/
OUTWEEKLY=$CINECA_SCRATCH/KdProduct/WEEKLYAVE/
SATMESH=SAT1km_mesh

mkdir -p $OUTMONTHLY
mkdir -p $TTOUTMONTHLY
mkdir -p $OUTWEEKLY


echo -----
echo mpirun -np 25 python aveSat_kd.py -i $CHECKED -o $OUTMONTHLY -m $SATMESH -t monthly

echo -----
echo mpirun -np 25 python aveSat_kd.py -i $CHECKED -o $OUTWEEKLY -m $SATMESH -t weekly_monday
