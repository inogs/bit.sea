export MASKFILE=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_09/wrkdir/MODEL/meshmask.nc

mkdir -p tmp_nc Fig_Hovmoeller_Stat
NCDIR=tmp_nc
OUTDIR=Fig_Hovmoeller_Stat
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR
