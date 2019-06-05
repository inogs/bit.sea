#DIRRST=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/DA_RST/
DIRRST=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/wrkdir/MODEL/DA__FREQ_1/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
#OUTCORR=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/MEANcorr/
OUTCORR=/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/DAdelta/MEANcorr/
mkdir -p $OUTCORR

echo python meancorrsub.py -i $DIRRST -m $MASKFILE -o $OUTCORR

