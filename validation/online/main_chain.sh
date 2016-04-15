OPA_RUNDATE=20160412

STARTTIME_a=$OPA_RUNDATE-10
STARTTIME_a=$OPA_RUNDATE-7

STARTTIME_f=$OPA_RUNDATE-7
END__TIME_f=$OPA_RUNDATE-4
MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc


ARCHIVE_DIR=/pico/home/usera07ogs/a07ogs00/OPA/V2C/archive

ONLINE_VALIDATION_DIR=/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data_02

mkdir -p $ONLINE_VALIDATION_DIR
python archive_extractor.py --type analysis -st ${STARTTIME_a} -et ${END__TIME_a}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}
python archive_extractor.py --type forecast -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}