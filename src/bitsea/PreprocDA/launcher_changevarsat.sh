

MESHMASK=/g100_scratch/userexternal/ateruzzi/PerCarolina/ADJUST_VARSAT/meshmask_v12c.nc
INDIR=/g100_scratch/userexternal/ateruzzi/PerCarolina/ADJUST_VARSAT/VARSAT_V9C
# LIMERRCOAST=0.1
LIMERRCOAST=nan
FIRSMONTH=6
LASTMONT=8
ADDERRSUMMER=0.005
SATMAP=Mesh24

OUTDIR=/g100_scratch/userexternal/ateruzzi/PerCarolina/ADJUST_VARSAT/VARSATadj_COAST${LIMERRCOAST}_SUMMER${ADDERRSUMMER}

mkdir -p $OUTDIR

echo python adjust_varsat -i $INDIR -o $OUTDIR -m $MESHMASK -s $FIRSMONTH -e $LASTMONT -v $ADDERRSUMMER -c $LIMERRCOAST -mm $SATMAP
