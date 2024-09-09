# from qualification_quid.sh in ../deliverables

#cfr time series surface chl
RUNref=DA_FLOAT_SAT/Summer/RUN_CR/
export MASKFILE=/pico/scratch/userexternal/ateruzzi/$RUNref/wrkdir/MODEL/meshmask.nc

DIRPKL=DA_FLOAT_SAT/Summer/
mkdir -p Fig4.2/Summer
python plot_timeseries_runs_open_Summer.py -i $DIRPKL -m $MASKFILE -o export_data_ScMYValidation_plan_open_sea.pkl -O Fig4.2/Summer/ -s 20150801 -e 20150907


#cfr profili runs
INDIR=$CINECA_SCRATCH/DA_FLOAT_SAT/Summer/
OUTDIR=PROFILEScomparison/Summer/
mkdir -p $OUTDIR
python compare_sim_profiles_smallsubs_Summer.py -i $INDIR -o $OUTDIR -s 20150801 -e 20150907 -m $MASKFILE

#cfr float runs
INPUTDIR=DA_FLOAT_SAT/Summer/
OUTDIR=TIMESER_FLOAT/Summer/
mkdir -p $OUTDIR
# Attention: use profiler_satfloats


python plot_timeseries_runs_floats_Summer.py -i $INPUTDIR -O $OUTDIR

python plot_timeseries_runs_floats_DAdates_Summer.py -i $INPUTDIR -O $OUTDIR

