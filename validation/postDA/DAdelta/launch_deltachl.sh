#DIRRST=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/DA_RST/
DIRRST=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/wrkdir/MODEL/DA__FREQ_1/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
#OUTNPY=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/MEANdelta/
OUTNPY=/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/DAdelta/MEANdelta/
mkdir -p $OUTNPY

echo python delta_chl.py -i $DIRRST -m $MASKFILE -o $OUTNPY

