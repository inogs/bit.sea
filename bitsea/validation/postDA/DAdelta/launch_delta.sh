DIRRST=/marconi_scratch/userexternal/ateruzzi/ELAB_RA_COAST/DAdelta/DA_RST/
MASKFILE=/marconi_work/OGS_dev_0/MULTIPLATFORM/meshmask.nc
OUTNPY=/marconi_scratch/userexternal/ateruzzi/ELAB_RA_COAST/DAdelta/MASSbal_txt/
mkdir -p $OUTNPY

echo python delta_vars.py -i $DIRRST -m $MASKFILE -v c -o $OUTNPY

