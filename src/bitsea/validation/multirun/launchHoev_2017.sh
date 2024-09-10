MASK=$CINECA_SCRATCH/MASKS24/meshmask.nc

RUN=HC_2017_assw
RUN=HC_2017_simdd
OUTHOEV=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/Hoevmoller/$RUN/
mkdir -p $OUTHOEV

#launch profiler
BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_HC2017/VALID_float/$RUN/PROFILATORE/
export MASKFILE=$MASK
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler

echo python profilerHC_2017.py # To be launched 1 time
# and without mpi4py (module unload mpi4py/....) 

# Time series of profiles quantities
STATSDIR=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/$RUN/STATS/
mkdir -p $STATSDIR
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $STATSDIR
echo python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $STATSDIR -o $OUTHOEV

# Hovmoller
echo python Hov_flots+model_vars.py -m $MASKFILE -o $OUTHOEV -i $STATSDIR

# Values of BIAS and RMSD al layers
OUTSTATS=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/STATS/$RUN/
mkdir -p $OUTSTATS
echo python biofloats_ms.py  -m $MASKFILE -o $STATSDIR/float_bias_rmse.nc
echo python biofloats_ms_plotter.py -i $STATSDIR/float_bias_rmse.nc -f $OUTSTATS -t $STATSDIR

# Statistics on profile quantities
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $STATSDIR
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $STATSDIR -o $OUTSTATS


# Comparison
BASESTATS=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/
OUTDIFF=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/Hoevmoller/DIFF/
mkdir -p $OUTDIFF
echo Hovdiff_models_vars.py -m $MASKFILE -i $BASESTATS -o $OUTDIFF


## Profiler for physics
BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_HC2017/VALID_float/$RUN/PROFILATORE_PHYS/
export MASKFILE=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_simdd/bin_phys/meshmask_INGV.nc
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler
echo python profiler_phys.py # To be launched 1 time

# Time series of profiles quantities
STATSDIR=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/$RUN/STATS/
mkdir -p $STATSDIR
# for physical variables
echo python SingleFloatPhys_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $STATSDIR

# Hovmoller
echo python HovPhys_flots+model_vars.py -m $MASKFILE -o $OUTHOEV -i $STATSDIR
echo python HovDens_flots+model_vars.py -m $MASKFILE -o $OUTHOEV -i $STATSDIR

# correlation and covariance
MASKBGC=$CINECA_SCRATCH/MASKS24/meshmask.nc
OUTFIGCORR=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/CorrCovn/
mkdir -p $OUTFIGCORR
echo python corr_profiles.py -p $MASKFILE -m $MASKBGC -o $STATSDIR
echo python plot_corr_profiles.py -i $BASESTATS -o $OUTFIGCORR

# Comparison
BASESTATS=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/
OUTDIFF=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/Hoevmoller/DIFF/
mkdir -p $OUTDIFF
echo python HovdiffPhys_models_vars.py -m $MASKFILE -i $BASESTATS -o $OUTDIFF
echo python HovdiffDens_models_vars.py -m $MASKFILE -i $BASESTATS -o $OUTDIFF
