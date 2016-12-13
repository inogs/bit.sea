# Comparison of time series of different runs

#RUN_LIST = [CR_COAST, DA_COAST_01]
#VAR_LIST = [P_i, ppn]
# Non so come passarli

export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc

python extract_satTimeSeries.py

SATDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/WEEKLY_V4/
OUTSAT=/pico/scratch/userexternal/ateruzzi/BITDOTSEA/bit.sea/validation/multirun/

MODDIR=$CINECA_SCRATCH/
TXTPATH=/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/
OUTBASE=$CINECA_SCRATCH/ELAB_DA_COAST/OUTPUTvalidation/CFRTIMESER/
OUTRUN=$OUTBASE/PKLFILES/
OUTDIR=$OUTBASE/FIGURES/

mkdir -p $OUTDIR
mkdir -p $OUTRUN

#python extract_satTimeSeries.py -o $OUTSAT -i $SATDIR -y 2013

python extract_runTimeSeries.py -o $OUTRUN -i $MODDIR -t $TXTPATH -r CR_COAST

python plot_cfrTimeSeries.py -o $OUTDIR -i $OUTRUN -s $OUTSAT

python stat_cfrTimeSeries.py -o $OUTBASE -i $OUTRUN 
