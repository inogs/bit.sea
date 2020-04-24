export ONLINE_REPO=/gpfs/work/OGS_prod_0/OPA/V6C/prod/inpdir/ONLINE/
export MASKFILE=/gpfs/work/OGS_prod_0/OPA/V6C/prod/wrkdir/analysis/2/MODEL/meshmask.nc
BASEDIR=/gpfs/work/OGS_prod_0/OPA/V6C/devel/inpdir/analysis/VALIDATION/FLOAT/PROFILATORE


mkdir -p tmp_nc Fig_Hovmoeller_Stat
NCDIR=tmp_nc
OUTDIR=Fig_Hovmoeller_Stat
RUNDATE=20200331
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR -d $RUNDATE

# python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -m $MASKFILE -b $BASEDIR -i $NCDIR -o $OUTDIR

echo python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR -d $RUNDATE
