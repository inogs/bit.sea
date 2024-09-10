# from qualification_quid.sh in ../deliverables

#cfr time series surface chl
RUNref=DA_FLOAT_SAT/RUN_CR_2015/ # per i run corti Winter Summer allruns_short.sh
export MASKFILE=/pico/scratch/userexternal/ateruzzi/$RUNref/wrkdir/MODEL/meshmask.nc

DIRPKL=DA_FLOAT_SAT/
mkdir -p Fig4.2/2015
python plot_timeseries_runs_open.py -i $DIRPKL -m $MASKFILE -o export_data_ScMYValidation_plan_open_sea.pkl -O Fig4.2/2015/ -s 20150101 -e 20151231


#cfr profili runs
INDIR=$CINECA_SCRATCH/DA_FLOAT_SAT/
OUTDIR=PROFILEScomparison/2015/
#INDIR=$CINECA_SCRATCH/DA_FLOAT_SAT/Winter/
#OUTDIR=PROFILEScomparison/Winter/
mkdir -p $OUTDIR


# LA PARTE SOTTO E' DA VEDERE 


#python compare_sim_profiles_smallsubs.py -i $INDIR -o $OUTDIR -s 20150101 -e 20151231 -m $MASKFILE
#python compare_sim_profiles_smallsubs_Winter.py -i $INDIR -o $OUTDIR -s 20150101 -e 20151231 -m $MASKFILE

#cfr float runs
#INPUTDIR=DA_FLOAT_SAT/Winter/
#OUTDIR=TIMESER_FLOAT/Winter/
#mkdir -p $OUTDIR
# Attention: use profiler_satfloats


#python plot_timeseries_runs_floats_Winter.py -i $INPUTDIR -O $OUTDIR

#python plot_timeseries_runs_floats_DAdates_Winter.py -i $INPUTDIR -O $OUTDIR

