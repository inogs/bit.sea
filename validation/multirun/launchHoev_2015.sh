MASK=$CINECA_SCRATCH/MASKS16/meshmask.nc

RUN=DA_FloatNut/RUN_chl_N3nN1pupd
RUN=DA_FloatNut/RUN_N3nN1p_chlfreq
RUN=DA_FloatNut/RUN_REF
RUN=DA_FloatNut/RUN_chl_1day
OUTHOEV=$CINECA_SCRATCH/ELAB_DAfloatNut/VALID_float/FIGURES/Hoevmoller/$RUN/
mkdir -p $OUTHOEV

#launch profiler
BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAfloatNut/VALID_float/$RUN/PROFILATORE/
export MASKFILE=$MASK
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler

echo python profiler_2015.py # To be launched 1 time
# and without mpi4py (module unload mpi4py/....) 

# Time series of profiles quantities
STATSDIR=$CINECA_SCRATCH/ELAB_DAfloatNut/VALID_float/$RUN/STATS/
mkdir -p $STATSDIR
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $STATSDIR
echo python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $STATSDIR -o $OUTHOEV

# Hovmoller
echo python Hov_flots+model_vars.py -m $MASKFILE -o $OUTHOEV -i $STATSDIR

# Values of BIAS and RMSD al layers
#OUTSTATS=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/STATS/$RUN/
OUTSTATS=$CINECA_SCRATCH/ELAB_DAfloatNut/VALID_float/FIGURES/STATS/$RUN/
mkdir -p $OUTSTATS
echo python biofloats_ms.py  -m $MASKFILE -o $STATSDIR/float_bias_rmse.nc
echo python biofloats_ms_plotter.py -i $STATSDIR/float_bias_rmse.nc -f $OUTSTATS -t $STATSDIR

module load scipy/1.1.0--python--2.7.12
# Statistics on profile quantities
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $STATSDIR
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $STATSDIR -o $OUTSTATS


# Comparison
BASESTATS=$CINECA_SCRATCH/ELAB_DAfloatNut/VALID_float/
DIRDIFF=$RUN
OUTDIFF=$CINECA_SCRATCH/ELAB_DAfloatNut/VALID_float/FIGURES/Hoevmoller/DIFF/$DIRDIFF
mkdir -p $OUTDIFF
echo Hovdiff_models_vars.py -m $MASKFILE -i $BASESTATS -o $OUTDIFF


