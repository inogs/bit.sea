INVARMOD=/gpfs/scratch/userexternal/ateruzzi/INTERP_VAR/VAR_MOD/VARMOD_16/
OUTVARMOD=/gpfs/scratch/userexternal/ateruzzi/INTERP_VAR/VAR_MOD/VARMOD_4/

MASK16=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
MASK4=/gpfs/scratch/userexternal/ateruzzi/MASKS4/meshmask.nc
FILEPREFIX=var2Dcoast.

mkdir -p $OUTVARMOD


echo python interp_var.py -i $INVARMOD -o $OUTVARMOD -m $MASK16 -n $MASK4 -p $FILEPREFIX



INVARSAT=/gpfs/scratch/userexternal/ateruzzi/INTERP_VAR/VAR_SAT/VARSAT_16/
OUTVARSAT=/gpfs/scratch/userexternal/ateruzzi/INTERP_VAR/VAR_SAT/VARSAT_4/

MASK16=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
MASK4=/gpfs/scratch/userexternal/ateruzzi/MASKS4/meshmask.nc
FILEPREFIX=var2D.

mkdir -p $OUTVARSAT


echo python interp_var.py -i $INVARSAT -o $OUTVARSAT -m $MASK16 -n $MASK4 -p $FILEPREFIX
