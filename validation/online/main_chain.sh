OPA_RUNDATE=20160426 # tuesday 

STARTTIME_a=$( date -d " $OPA_RUNDATE -10 days " +%Y%m%d ) 
END__TIME_a=$( date -d " $OPA_RUNDATE  -8 days " +%Y%m%d )

STARTTIME_s=$( date -d " $OPA_RUNDATE -7 days " +%Y%m%d )

STARTTIME_f=$( date -d " $OPA_RUNDATE  -7 days " +%Y%m%d )
END__TIME_f=$( date -d " $OPA_RUNDATE  -4 days " +%Y%m%d )

MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc

WRKDIR=/pico/home/usera07ogs/a07ogs00/OPA/V2C/wrkdir/5
ARCHIVE_DIR=/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/archive

ONLINE_VALIDATION_DIR=/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data_V2
TMP_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/TMP
PROFILERDIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/PROFILATORE
IMG_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/IMG/

TMP_DIR__ACT=${ONLINE_VALIDATION_DIR}/ACTUAL/TMP
PROFILERDIRA=${ONLINE_VALIDATION_DIR}/ACTUAL/PROFILATORE
IMG_DIR__ACT=${ONLINE_VALIDATION_DIR}/ACTUAL/IMG/

mkdir -p $ONLINE_VALIDATION_DIR/PREVIOUS
mkdir -p $ONLINE_VALIDATION_DIR/ACTUAL

if [ 1 == 0 ] ; then
python archive_extractor.py --type analysis -st ${STARTTIME_a} -et ${END__TIME_a}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_s} -et ${STARTTIME_s}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type analysis -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/ACTUAL
# ho ottenuto una serie continua dall'archivio, ma analisi e forecasts della precedente
#step 2 -- re-do Forcings generation

./forcings_gen.sh -d ${ONLINE_VALIDATION_DIR}/PREVIOUS
./forcings_gen.sh -d ${ONLINE_VALIDATION_DIR}/ACTUAL
# step 3 -- aggregate variables to have chl and all variables (physical and biological) on the same ave files

python var_aggregator.py --biodir ${ONLINE_VALIDATION_DIR}/PREVIOUS/output_bio --physdir ${ONLINE_VALIDATION_DIR}/PREVIOUS/FORCINGS -t $TMP_DIR_PREV

python var_aggregator.py --biodir ${ONLINE_VALIDATION_DIR}/ACTUAL/output_bio --physdir ${ONLINE_VALIDATION_DIR}/ACTUAL/FORCINGS -t $TMP_DIR__ACT



SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/WEEKLY/
SAT_DAILY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/
PREVIOUS_TUE_RUNDIR=${STARTTIME_s}
f0_name=${ONLINE_VALIDATION_DIR}/Validation_f0_${PREVIOUS_TUE_RUNDIR}_on_weekly_Sat.${STARTTIME_s}.nc
python SatValidation.py -d ${STARTTIME_s} -f $TMP_DIR_PREV -s ${SAT_WEEKLY_DIR} -o $f0_name  # MISFIT   

for fc_day in 1 2; do
   DAY=$(date -d "$STARTTIME_s + $fc_day days"  +%Y%m%d )
   f_name=${ONLINE_VALIDATION_DIR}/Validation_f${fc_day}_${PREVIOUS_TUE_RUNDIR}_on_daily_Sat.${DAY}.nc
   python SatValidation.py  -d $DAY -f $TMP_DIR_PREV -s ${SAT_DAILY_DIR} -o ${f_name}     # ERROR  
done

# Questi devono essere archiviati in $ARCHIVE_DIR/POSTPROC/AVE_FREQ_1/validation/Sat


a0_name=${ONLINE_VALIDATION_DIR}/Validation_a0_${OPA_RUNDATE}_on_weekly_Sat.${STARTTIME_s}.nc
python SatValidation.py -d ${STARTTIME_s} -f $TMP_DIR__ACT -s ${SAT_WEEKLY_DIR} -o $a0_name  # MISFIT   

for ac_day in 1 2; do
   DAY=$(date -d "$STARTTIME_s + $ac_day days"  +%Y%m%d )
   a_name=${ONLINE_VALIDATION_DIR}/Validation_a${ac_day}_${OPA_RUNDATE}_on_daily_Sat.${DAY}.nc
   python SatValidation.py  -d $DAY -f $TMP_DIR__ACT -s ${SAT_DAILY_DIR} -o ${a_name}     # ERROR  
done


fi

#wrkdir=/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/wrkdir/5

# python adding_vars.py -v vosaline -i $WRKDIR/MODEL/FORCINGS -ip T* --input_prefix T -o $wrkdir/POSTPROC/AVE_FREQ_1/TMP -op ave*nc --output_prefix ave. -st 20160617 -et 20160619 -m $MASKFILE
# python adding_vars.py -v votemper -i $WRKDIR/MODEL/FORCINGS -ip T* --input_prefix T -o $wrkdir/POSTPROC/AVE_FREQ_1/TMP -op ave*nc --output_prefix ave. -st 20160617 -et 20160619 -m $MASKFILE
#python var_aggregator.py --biodir $WRKDIR/POSTPROC/AVE_FREQ_1/TMP --physdir $WRKDIR/MODEL/FORCINGS -t $TMP_DIR__ACT

exit 0


# Matchup analysis
python float_extractor.py -st ${STARTTIME_a} -et ${END__TIME_a} -i $TMP_DIR_PREV -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
#Matchup forecasts
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i $TMP_DIR_PREV -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
# Matchup actual
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i $TMP_DIR__ACT -b $PROFILERDIRA  -o $IMG_DIR__ACT -m $MASKFILE


#${STARTTIME_a} -et ${END__TIME_a} confronto Float-Analysis previous
#${STARTTIME_f} -et ${END__TIME_f} confronto Float con {analisi attuale, IMG_DIR__ACT}, {forecast, IMG_DIR_PREV}
python profile_plotter3.py -f $IMG_DIR_PREV -a $IMG_DIR__ACT # scorre tutti i files e genera le immagini a 3


Bisogna fare tutti i link !!!!!
python nrt3.py -o $OUTFILE -b $BASEDIR -m $MASKFILE -i $INPUTDIR/TMP -d $DATE

ARCHIVEDIR=/pico/home/usera07ogs/a07ogs00/OPA/V2C/archive
ARCHIVE_PREV=/pico/scratch/userexternal/gbolzon0/NRT/V4/NRT3_outputs/


VALID_WRKDIR=/pico/scratch/userexternal/gbolzon0/NRT/V2C/NRT3_outputs/
OUTDIR=IMG
mkdir -p $OUTDIR

python nrt3_plotter.py -a $ARCHIVEDIR -p $ARCHIVE_PREV -v $VALID_WRKDIR -o IMG/ -d $STARTTIME_s