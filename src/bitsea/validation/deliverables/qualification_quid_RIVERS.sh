# Descriptor of CMEMS-Med-biogeochemistry-ScQP-v3.1.docx

# QUID FOR RIVERS VALIDATION:
# SECTION 4:
. ../online/profile.inc

export MASKFILE=/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v3/wrkdir/MODEL/meshmask.nc 
export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/

         INPUTDIR=/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v3/wrkdir/MODEL/AVE_FREQ_2/
   INPUT_AGGR_DIR=/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v3/wrkdir/MODEL/AVE_FREQ_2/
STAT_PROFILES_DIR=/g100_scratch/userexternal/gbolzon0/V10C/run4.19/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/

SAT_KD_WEEKLY_DIR=/g100_scratch/userexternal/lfeudale/KD490/SAT/KD490/DT/SIGMA_2.0_KDmin0.022/WEEKLY_4_24/
SAT_CHLWEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C/SAT/CHL/DT/WEEKLY_4_24/
SAT_CHLPFTsWEEKLY_DIR=/g100_scratch/userexternal/gbolzon0/V10C/SAT/CHL/DT/WEEKLY_4_24
KD_MODEL_DIR=/g100_scratch/userexternal/gbolzon0/V10C/run4.19/wrkdir/POSTPROC/output/AVE_FREQ_3/KD_WEEKLY/
BASEDIR=/g100_scratch/userexternal/lfeudale/validation/V10C/run4.19/PROFILATORE/


# CREATE Directories for SAT:
#opa_prex_or_die "mkdir -p Fig4.2 Fig4.3/offshore Fig4.3/coast table4.1 table4.2"
opa_prex_or_die "mkdir -p ./SAT_rivers/Timeseries/ ./SAT_rivers/BIAS_RMSD/"


#for VAR in   kd490  P_l  P1l    P2l    P3l    P4l  ; do
for VAR in   P_l  ; do
    RIVERS_PKL={$VAR}_rivers.pkl

    MODELDIR=$INPUTDIR
    SAT_DIR=$SAT_CHLWEEKLY_DIR

    if [ ${VAR} == 'P_l' ] ; then 
        MODELDIR=$INPUT_AGGR_DIR #; fi
        SAT_DIR=$SAT_CHLWEEKLY_DIR
     fi

    if [ $VAR == 'kd490' ] ; then 
        MODELDIR=$KD_MODEL_DIR
        SAT_DIR=$SAT_KD_WEEKLY_DIR
     fi
 
    if [$VAR == 'P1l'] | [$VAR == 'P2l'] | [$VAR == 'P3l'] | [$VAR == 'P4l']; then 
    # Plot PFTs comparison with HPLC dataset:
       INDIR_HPLC_CLIM=/g100_scratch/userexternal/lfeudale/validation/V10C/run4.19/bit.sea/validation/deliverables/
       opa_prex_or_die "python plot_HPLC_clim.py -i $INDIR_HPLC_CLIM -p $VAR "
     fi

    LAYER=10
   # if [ $VAR == 'kd490' ] ; then LAYER=0 ; fi
   
   
    opa_prex_or_die "python ScMYvalidation_plan_RIVERS.py -v $VAR -s $SAT_DIR -i $MODELDIR -m $MASKFILE -c everywhere -l $LAYER  -o $RIVERS_PKL --datestart 20190101 --dateend 20200101"

    opa_prex_or_die "python plot_timeseries_STD_RIVERS.py -v $VAR -i $RIVERS_PKL -o ./SAT_rivers/Timeseries/ "
    opa_prex_or_die "python plot_timeseries_RMS_CORR_RIVERS.py -v $VAR -i $RIVERS_PKL -o ./SAT_rivers/BIAS_RMSD/"   # table4.2

done

exit 0
