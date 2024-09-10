GSSOBS_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/


OUT_STATS=/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/CFR_CLIMA/PKL/
mkdir -p $OUT_STATS
CLIM_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SUBBASIN_STATISTICS/
SAT_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/TIMESER_SATMONTHLY_RACOAST/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
NAMEFILEOUT=rea16composed

echo python climasat_timeRA.py -c $CLIM_DIR -m $MASKFILE  -s $SAT_DIR -o $OUT_STATS -f $NAMEFILEOUT


echo -----------

OUT_FIG=/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/CFR_CLIMA/FIG/
mkdir -p $OUT_FIG

echo python plot_climasat_timeRA.py  -s $OUT_STATS -o $OUT_FIG

echo -----------


