DIRRST=/marconi_scratch/userexternal/ateruzzi/ELAB_RA_COAST/DAdelta/DA_RST/
MASKFILE=/marconi_work/OGS_dev_0/MULTIPLATFORM/meshmask.nc
OUTMISF=/marconi_scratch/userexternal/ateruzzi/ELAB_RA_COAST/DAdelta/MEANmisf/
mkdir -p $OUTMISF

echo python meanmisf.py -i $DIRRST -m $MASKFILE -o $OUTMISF

