# from qualification_quid.sh in ../deliverables

#cfr time series surface chl
RUNref=DA_FLOAT_SAT/Winter/RUN_SAT_01/
export MASKFILE=/pico/scratch/userexternal/ateruzzi/$RUNref/wrkdir/MODEL/meshmask.nc

DIRPKL=DA_FLOAT_SAT/Winter/
mkdir -p Fig4.2/Winter
python plot_timeseries_runs_open.py -i $DIRPKL -m $MASKFILE -o export_data_ScMYValidation_plan_open_sea.pkl -O Fig4.2/Winter/

#cfr profili runs
INDIR=$CINECA_SCRATCH/DA_FLOAT_SAT/Winter/
OUTDIR=PROFILEScomparison/Winter/
mkdir -p $OUTDIR
python compare_sim_profiles_smallsubs.py -i $INDIR -o $OUTDIR -s 20150101 -e 20150209 -m $MASKFILE

#cfr float runs
INPUTDIR=DA_FLOAT_SAT/Winter/
OUTDIR=TIMESER_FLOAT/Winter/
mkdir -p $OUTDIR
# Attention: use profiler_satfloats


python plot_timeseries_runs_floats.py -i $INPUTDIR -O $OUTDIR

python plot_timeseries_runs_floats_DAdates.py -i $INPUTDIR -O $OUTDIR

