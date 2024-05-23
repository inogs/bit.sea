#! /bin/bash

. ../../Sat/profile.inc

export ONLINE_REPO=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/DA
export PROFILES_SOURCE=ppcon

DA_TIMES_SAT=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/run4.0/wrkdir/MODEL/daTimes_sat
DAFREQ=1
DA_HH=13

DATIMES_WRKDIR=./
OUTDIR=DA_NML_NEW
OUTFILE=${DATIMES_WRKDIR}/daTimes


mkdir -p $OUTDIR $DATIMES_WRKDIR

for VARDA in N3n P_l O2o ; do
   my_prex_or_die "python profile_dates_DAfreq.py -f $DAFREQ -t $DA_HH -s 20190101 -e 20200101 -v ${VARDA} -o ${DATIMES_WRKDIR}/daTimes_float_${VARDA} "
done


# Edit varlist in merge_daTimes.py
my_prex_or_die "python merge_daTimes.py -s $DA_TIMES_SAT -d $DATIMES_WRKDIR -o ${OUTFILE}"


my_prex_or_die "python gen_3dvar_satfloat_namelists.py -o $OUTDIR -t satfloat.template.nml -r $OUTFILE -s $DA_TIMES_SAT -d $DATIMES_WRKDIR"
