#! /bin/bash

export ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/ONLINE_V5C/
DA_TIMES_SAT=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_04/wrkdir/MODEL/daTimes_sat
DAFREQ=1
DA_HH=13

DATIMES_WRKDIR=./
OUTDIR=DA_NML_NEW
OUTFILE=${DATIMES_WRKDIR}/daTimes


mkdir -p $OUTDIR $DATIMES_WRKDIR

for VARDA in N3n P_l; do
echo python profile_dates_DAfreq.py -f $DAFREQ -t $DA_HH -s 20180101 -e 20190101 -v ${VARDA} -o ${DATIMES_WRKDIR}/daTimes_float_${VARDA}
done


# Edit varlist in merge_daTimes.py
echo python merge_daTimes.py -s $DA_TIMES_SAT -d $DATIMES_WRKDIR -o ${OUTFILE}


python gen_3dvar_satfloat_namelists.py -o $OUTDIR -t satfloat.template.nml -r $OUTFILE -s $DA_TIMES_SAT -d $DATIMES_WRKDIR