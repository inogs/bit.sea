MASK=$CINECA_SCRATCH/MASKS16corrected/meshmask.nc

RUN=DA_Float/RUN_FLOAT_chl12
RUN=DA_Float/RUN_FLOAT_chl_nupd
RUN=DA_Float/RUN_FLOAT_chl_n
RUN=DA_Float/RUN_FLOAT_chlnit
RUN=DA_Float/RUN_FLOAT_nit
RUN=DA_Float/RUN_FLOAT_chlnit_std2d
OUTHOEV=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/HoevmollerINC/$RUN/
mkdir -p $OUTHOEV

#launch profiler
export MASKFILE=$MASK
export ONLINE_REPO=/gpfs/scratch/userexternal/ateruzzi/REPO_LOV/ONLINE/ # THE SAME OF DA
BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/$RUN/PROFILATORE_RSTaft/
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler

echo python profilerRSTaft_2015.py # To be launched 1 time
# and without mpi4py (module unload mpi4py/....) 

BASEDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/$RUN/PROFILATORE_RSTbef/
mkdir -p $BASEDIR # Must be consistent with BASEDIRE in profiler
echo python profilerRSTbef_2015.py # To be launched 1 time
# and without mpi4py (module unload mpi4py/....) 


# Hovmoeller RST
OUTHOVRST=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/HoevmollerRST/$RUN/
mkdir -p $OUTHOVRST
echo python Hov_flots+model_varsRSTbef.py -m $MASKFILE -o $OUTHOVRST
echo python Hov_flots+model_varsRSTaft.py -m $MASKFILE -o $OUTHOVRST

# Hovmoller increments
echo python Hov_DAincrements.py -m $MASKFILE -o $OUTHOEV



OUTHOEVINN=$CINECA_SCRATCH/ELAB_DAFloat/VALID_float/FIGURES/HoevmollerINNOV/$RUN/
mkdir -p $OUTHOEVINN
FLOATPREP=$CINECA_SCRATCH/$RUN/wrkdir/float_preproc/
# Hovmoller innnovation
echo python Hov_DAinnovation.py -m $MASKFILE -o $OUTHOEVINN

