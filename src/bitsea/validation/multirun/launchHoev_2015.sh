MASK=$CINECA_SCRATCH/MASKS16corrected/meshmask.nc

RUN=DA_Float/RUN_SAT_FLOAT_chl_n
OUTHOEV=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/Hoevmoller/$RUN/
OUTHOEVMODEL=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/HoevmollerModel/$RUN/
mkdir -p $OUTHOEV
mkdir -p $OUTHOEVMODEL

#launch profiler
BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/$RUN/PROFILATORE/
export MASKFILE=$MASK
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler

echo python profiler_2015.py # To be launched 1 time
# and without mpi4py (module unload mpi4py/....) 

# Time series of profiles quantities
STATSDIR=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/$RUN/STATS/
mkdir -p $STATSDIR
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $STATSDIR
echo python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $STATSDIR -o $OUTHOEV

# Hovmoller
echo python Hov_flots+model_vars.py -m $MASKFILE -o $OUTHOEV -i $STATSDIR
# Hovmoeller of model variables at float positions
echo python Hov_onlymodel_vars.py -m $MASKFILE -o $OUTHOEVMODEL


# Values of BIAS and RMSD al layers
#OUTSTATS=$CINECA_SCRATCH/ELAB_HC2017/VALID_float/FIGURES/STATS/$RUN/
OUTSTATS=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/STATS/$RUN/
mkdir -p $OUTSTATS
echo python biofloats_ms.py  -m $MASKFILE -o $STATSDIR/float_bias_rmse.nc
echo python biofloats_ms_plotter.py -i $STATSDIR/float_bias_rmse.nc -f $OUTSTATS -t $STATSDIR

module load scipy/1.1.0--python--2.7.12
# Statistics on profile quantities
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $STATSDIR
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $STATSDIR -o $OUTSTATS


# Comparison
BASESTATS=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/
DIRDIFF=$RUN
OUTDIFF=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/Hoevmoller/DIFF/$DIRDIFF
mkdir -p $OUTDIFF
echo Hovdiff_models_vars.py -m $MASKFILE -i $BASESTATS -o $OUTDIFF

#Hov density
BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/$RUN/PROFILATORE_PHYS/
OUTHPHYS=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/$RUN/FIGURES/HovmoellerPHY/
mkdir -p $OUTHPHYS
export MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16_INGV/meshmask_1672.nc
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler
echo python profiler_phys.py # To be launched 1 time

echo python HovDens_flots+model_vars.py -m $MASKFILE -o $OUTHPHYS


#Hov N with slope
MASKFILE=$CINECA_SCRATCH/MASKS16corrected/meshmask.nc
OUTSLOPE=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/$RUN/FIGURES/HovmoellerNslope/
mkdir -p $OUTSLOPE
echo python Hov_flots+model_slopeN.py -m $MASKFILE -o $OUTSLOPE 
