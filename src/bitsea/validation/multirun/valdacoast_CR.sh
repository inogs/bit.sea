# Validatio of coastal assimilation runs

export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/WEEKLY_V4/

export START_DATE=20130101
export END_DATE=20140101

RUN=CR_COAST

MOD_DAILY_DIR=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/
MOD_MONTHLY_DIR=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/

OUTDIR=/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/OUTPUTvalidation/$RUN
mkdir $OUTDIR

# Time series
mkdir -p $OUTDIR/TIMESER/
python ScMYvalidation_plan.py -o $OUTDIR/export_data_ScMYValidation_plan.pkl -s $SAT_WEEKLY_DIR -m $MOD_DAILY_DIR -l 0.

python plot_timeseries.py -i $OUTDIR/export_data_ScMYValidation_plan.pkl -o $OUTDIR/TIMESER/ -r $RUN


# RMS and BIAS time series
mkdir -p $OUTDIR/BIAS_RMS/
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan.pkl -o $OUTDIR/BIAS_RMS/ -r $RUN


# Media annuale
# ! year defined with START_DATE and END_DATE 
mkdir -p $OUTDIR/MEAN_MAPS

python averager_and_plot_map.py -i $MOD_MONTHLY_DIR  -o $OUTDIR/MEAN_MAPS/ -v P_i -t mean -l multilayer # Layers hardcoded
python sat_ave_and_plot.py      -i $SAT_WEEKLY_DIR  -o $OUTDIR/MEAN_MAPS/

python averager_and_plot_map.py -i $MOD_MONTHLY_DIR  -o $OUTDIR/MEAN_MAPS/ -v N1p -t mean -l multilayer # Layers hardcoded
python averager_and_plot_map.py -i $MOD_MONTHLY_DIR  -o $OUTDIR/MEAN_MAPS/ -v N3n -t mean -l multilayer # Layers hardcoded

# PPN
mkdir -p $OUTDIR/PPN
python averager_and_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/PPN/ -v ppn -t integral -l 0-200layer


# mappe mensili chl
mkdir -p $OUTDIR/MONTHLY_MAPS

python monthly_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/MONTHLY_MAPS/  -v P_i -t mean # Layers hardcoded
python sat_monthly_plot_map.py -i $SAT_WEEKLY_DIR -o $OUTDIR/MONTHLY_MAPS/

