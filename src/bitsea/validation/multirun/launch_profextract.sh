MASK=$CINECA_SCRATCH/MASKS16/meshmask.nc

RUN=DA_FloatNut/RUN_chl_1day
RUN=DA_FloatNut/RUN_N3nN1p_chlfreq
RUN=DA_FloatNut/RUN_chl_N3nN1pupd


# profiler should have been already executed on the runs
# listed in the profiler

BASESTATS=$CINECA_SCRATCH/ELAB_DAfloatNut/VALID_float/
DIRPROF=$RUN
OUTPROF=$CINECA_SCRATCH/ELAB_DAfloatNut/PROFILESextr/$DIRPROF
WMO=6901771

mkdir -p $OUTPROF

echo python profscomparison.py -m $MASK -i $BASESTATS -o $OUTPROF -w $WMO

