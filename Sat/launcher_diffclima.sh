INCLIMA=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SatClimatology.nc
INSAT=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C/SAT/MULTISENSOR/1Km/NRT/DAILY/CHECKED/
OUTDIR=/gpfs/scratch/userexternal/ateruzzi/SAT_DIFFCLIMA_2019/
OUTDIR=/gpfs/scratch/userexternal/ateruzzi/SAT_DIFFCLIMA/
STARTTIME=20190801
ENDTIME=20190831
STARTTIME=20200816
ENDTIME=20200831

mkdir -p $OUTDIR/DIFF
mkdir -p $OUTDIR/DIFF_N

echo python sat_clim_diff.py -i $INSAT -c $INCLIMA -o $OUTDIR -s $STARTTIME -e $ENDTIME -m SAT1km_mesh
