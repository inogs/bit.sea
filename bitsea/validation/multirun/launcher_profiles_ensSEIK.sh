#! /bin/bash

MASKFILE=/g100_scratch/userexternal/ateruzzi/PerSimone/EVAL_SEIK/POSTPROC/meshmask.nc
RUN=RST0101_ensCHL_3rd
INDIR=/g100_scratch/userexternal/ateruzzi/PerSimone/EVAL_SEIK/ENS_IC/out${RUN}/
OUTDIR=/g100_scratch/userexternal/ateruzzi/PerSimone/EVAL_SEIK/ENS_IC/PROFILES_ENS/${RUN}/

mkdir -p $OUTDIR

echo python profiles_plotter_ensSEIK.py -i $INDIR -o $OUTDIR -m $MASKFILE
