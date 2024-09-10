#DIRRST=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/DA_RST/
DIRRST=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/wrkdir/MODEL/DA__FREQ_1/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
#OUTMISF=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/MEANmisf/
OUTMISF=/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/DAdelta/MEANmisf/
mkdir -p $OUTMISF

echo python meanmisfsub.py -i $DIRRST -m $MASKFILE -o $OUTMISF

