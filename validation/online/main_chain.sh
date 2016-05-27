OPA_RUNDATE=20160525

STARTTIME_a=$( date -d " $OPA_RUNDATE -10 days " +%Y%m%d ) 
END__TIME_a=$( date -d " $OPA_RUNDATE  -7 days " +%Y%m%d )

STARTTIME_f=$( date -d " $OPA_RUNDATE  -7 days " +%Y%m%d )
END__TIME_f=$( date -d " $OPA_RUNDATE  -7 days " +%Y%m%d )
MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc


ARCHIVE_DIR=/pico/home/usera07ogs/a07ogs00/OPA/V2C/archive

ONLINE_VALIDATION_DIR=/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data_02
TMP_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/TMP
PROFILERDIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/PROFILATORE
IMG_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/IMG/

TMP_DIR__ACT=${ONLINE_VALIDATION_DIR}/ACTUAL/TMP
PROFILERDIRA=${ONLINE_VALIDATION_DIR}/ACTUAL/PROFILATORE
IMG_DIR__ACT=${ONLINE_VALIDATION_DIR}/ACTUAL/IMG/

mkdir -p $ONLINE_VALIDATION_DIR/PREVIOUS

python archive_extractor.py --type analysis -st ${STARTTIME_a} -et ${END__TIME_a}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS

# ho ottenuto una serie continua dall'archivio, ma analisi e forecasts della precedente
#step 2 -- re-do Forcings generation
./forcings_gen.sh -d ${ONLINE_VALIDATION_DIR}/PREVIOUS

# step 3 -- aggregate variables to have chl and all variables (physical and biological) on the same ave files

python var_aggregator.py --biodir ${ONLINE_VALIDATION_DIR}/output_bio --physdir ${ONLINE_VALIDATION_DIR}/FORCINGS -t $TMP_DIR_PREVIOUS


python var_aggregator.py --biodir wrkdir/2/POSTPROC/AVE_FREQ_1/TMP --physdir wrkdir/2/MODEL/FORCINGS -t $TMP_DIR_ACTUAL

# Matchup analysis
python float_extractor.py -st ${STARTTIME_a} -et ${END__TIME_a} -i $TMP_DIR_PREV -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
#Matchup forecasts
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i $TMP_DIR_PREV -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
# Matchup actual
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i $TMP_DIR__ACT -b $PROFILERDIRA  -o $IMG_DIR__ACT -m $MASKFILE


#${STARTTIME_a} -et ${END__TIME_a} confronto Float-Analysis previous
#${STARTTIME_f} -et ${END__TIME_f} confronto Float con {analisi attuale, IMG_DIR__ACT}, {forecast, IMG_DIR_PREV}
python profile_plotter3.py -f $IMG_DIR_PREV -a $IMG_DIR__ACT # scorre tutti i files e genera le immagini a 3

