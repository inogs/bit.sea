CHECKED_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST/2018/SAT_CORRECTED/
MONTHLY_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST/MULTI_MONTHLY2018/
DIR_DATE=/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST/DATES_MONTHLY2018/
mkdir -p $MONTHLY_DIR
mkdir -p $DIR_DATE

echo mpirun -np 10 python aveSat.py -i $CHECKED_DIR -o $MONTHLY_DIR -m V4mesh -t monthly -d $DIR_DATE

