#! /bin/bash

#SBATCH --job-name=daTimes
#SBATCH -N1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:30:00
#SBATCH --mem=100gb
#SBATCH --account=OGS_devC
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo

source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh
unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false
. ../../Sat/profile.inc

export ONLINE_REPO=/g100_work/OGS_devC/V11C/TRANSITION/ONLINE
export PROFILES_SOURCE=ppcon

DA_TIMES_SAT=/g100_work/OGS_devC/V11C/TRANSITION/wrkdir/MODEL/daTimes_sat
DAFREQ=1
DA_HH=13

DATIMES_WRKDIR=./
OUTDIR=DA_NML_NEW
OUTFILE=${DATIMES_WRKDIR}/daTimes


mkdir -p $OUTDIR $DATIMES_WRKDIR

for VARDA in N3n P_l O2o ; do
   my_prex_or_die "python profile_dates_DAfreq.py -f $DAFREQ -t $DA_HH -s 20210101 -e 20250101 -v ${VARDA} -o ${DATIMES_WRKDIR}/daTimes_float_${VARDA} "
done


# Edit varlist in merge_daTimes.py
my_prex_or_die "python merge_daTimes.py -s $DA_TIMES_SAT -d $DATIMES_WRKDIR -o ${OUTFILE}"


my_prex_or_die "python gen_3dvar_satfloat_namelists.py -o $OUTDIR -t satfloat.template.nml -r $OUTFILE -s $DA_TIMES_SAT -d $DATIMES_WRKDIR"
