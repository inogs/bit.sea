export MASKFILE=/gpfs/work/OGS_prod_0/OPA/V5C/devel/wrkdir/2/MODEL/meshmask.nc
BASEDIR=/gpfs/work/OGS_prod_0/OPA/V5C/prod/inpdir/VALIDATION/FLOAT/PROFILATORE


mkdir -p tmp_nc Fig_Hovmoeller_Stat
NCDIR=tmp_nc
OUTDIR=Fig_Hovmoeller_Stat
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
echo python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR
python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR

