#! /bin/bash

#PBS -N MAPS
#PBS -l walltime=0:30:00
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=100GB
#PBS -q route
#PBS -A OGS_dev_0

# cd $PBS_O_WORKDIR
RUNDATE=20190506
MASKFILE=/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/meshmask.nc

V4DIR=/gpfs/scratch/userexternal/gbolzon0/TRANSITION_24/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES
V5DIR=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_09/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES

TARFILE=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_09/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES.tar
BASEDIR=/gpfs/scratch/userexternal/gbolzon0/TRANSITION_24/wrkdir/POSTPROC/
  V5DIR=$BASEDIR/STAT_PROFILES
cd $BASEDIR
tar -xf $TARFILE -C $BASEDIR
python archive_aveScan_extractor.py -d $RUNDATE -a /gpfs/work/OGS_prod_0/OPA/V5C/prod/archive/ -o $V5DIR
mpirun -np 15 python compact_profiles -i $V5DIR -o $V5DIR



IMGDIR=/gpfs/scratch/userexternal/gbolzon0/TRANSITION_24/wrkdir/POSTPROC/IMG
mkdir -p $IMGDIR
mpirun -np 12 python profiles_plotter.py -V4 $V4DIR -V5 $V5DIR -m $MASKFILE -o $IMGDIR -d $RUNDATE -l ../deliverables/Plotlist_bio.xml




ARCHIVE_DIR=/gpfs/work/OGS_prod_0/OPA/V5C/devel/archive
  ZIPPED_DIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/MONTHLY_PRODUCTS/DAILY_gz
UNZIPPED_DIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/MONTHLY_PRODUCTS/DAILY
 MONTHLY_DIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/MONTHLY_PRODUCTS/MONTHLY_AVE
MONTHLY_PROD=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/MONTHLY_PRODUCTS/MONTHLY_PROD

mkdir -p $ZIPPED_DIR $UNZIPPED_DIR $MONTHLY_DIR $MONTHLY_PROD
POSTPROC_DIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/postproc/
python archive_monthly_extractor.py -a $ARCHIVE_DIR -o $ZIPPED_DIR
mpirun -np 20 python $POSTPROC_DIR/archive/uncompress.py -i $ZIPPED_DIR -o $UNZIPPED_DIR -l *gz


mpirun -np 10 python monthly_averager.py -i $UNZIPPED_DIR -o $MONTHLY_DIR -m $MASKFILE


basename $( ls $MONTHLY_DIR/ave*N1p.nc ) |  cut -c 5-12 > timelist.txt

RUNDATE=20190410
mpirun -np 1 python ${POSTPROC_DIR}/prodotti/prodotti_copernicus.py -i $MONTHLY_DIR -o $MONTHLY_PROD -t timelist.txt -d an -m $MASKFILE --tr monthly -b $RUNDATE







if [ 1 == 0 ]; then
module purge
module load profile/advanced
module load autoload
module load intel/pe-xe-2017--binary
module load netcdf/4.4.1--intel--pe-xe-2017--binary
module load python/2.7.12 scipy/0.18.1--python--2.7.12

source /marconi_work/OGS_dev_0/COPERNICUS/py_env_2.7.12/bin/activate
PYTHONPATH=$PYTHONPATH:/marconi/home/usera07ogs/a07ogs01/MAPPE/bit.sea
module load mpi4py/2.0.0--python--2.7.12

fi
#######################################
. profile.inc
#######################################

OPA_RUNDATE=20171128 # tuesday 

STARTTIME_a=$( date -d " $OPA_RUNDATE -10 days " +%Y%m%d ) 
END__TIME_a=$( date -d " $OPA_RUNDATE  -8 days " +%Y%m%d )

STARTTIME_s=$( date -d " $OPA_RUNDATE -7 days " +%Y%m%d )

STARTTIME_f=$( date -d " $OPA_RUNDATE  -7 days " +%Y%m%d )
END__TIME_f=$( date -d " $OPA_RUNDATE  -4 days " +%Y%m%d )

export MASKFILE=/marconi/home/usera07ogs/a07ogs00/OPA/V3C/etc/static-data/MED24_125/meshmask.nc
      INGV_MASK=/marconi/home/usera07ogs/a07ogs00/OPA/V3C/etc/static-data/MED24_141/meshmask.nc
WRKDIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/wrkdir/2
ARCHIVE_DIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/archive

ONLINE_VALIDATION_DIR=/marconi_scratch/usera07ogs/a07ogs01/online_validation_data/
ONLINE_VALIDATION_DIR_INPUT=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/wrkdir/2/POSTPROC/AVE_FREQ_1/online_validation/

if [ 1 == 0 ]; then
mkdir -p ${ONLINE_VALIDATION_DIR}/ANALYSIS_PROFILES_1week/PROFILATORE/PROFILES ${ONLINE_VALIDATION_DIR}/biofloat_ms
#cp Archive_biofloats_ms_validation_V2C.tar
for I in `ls ${ONLINE_VALIDATION_DIR_INPUT}/PREVIOUS/PROFILATORE/PROFILES/*nc ` ; do
  ln -fs $I ${ONLINE_VALIDATION_DIR}/ANALYSIS_PROFILES_1week/PROFILATORE/PROFILES
done
for I in `ls ${ONLINE_VALIDATION_DIR_INPUT}/ACTUAL/PROFILATORE/PROFILES/*nc  ` ; do
  ln -fs $I ${ONLINE_VALIDATION_DIR}/ANALYSIS_PROFILES_1week/PROFILATORE/PROFILES
done

NRT_biofloat_out_name=${ONLINE_VALIDATION_DIR}/ANALYSIS_PROFILES_1week/NRT_validation_${OPA_RUNDATE}_on_week_floats.${STARTTIME_s}.nc

python biofloats_ms.py -d $STARTTIME_s -m $MASKFILE -b ${ONLINE_VALIDATION_DIR}/ANALYSIS_PROFILES_1week/PROFILATORE -o ${NRT_biofloat_out_name}
fi

#-----------------
OLD_ARCHIVE=Archive_biofloats_ms_validation_V2.tar
mkdir -p ${ONLINE_VALIDATION_DIR}/biofloat_ms/V2
#cp ${OPA_ETCDIR}/static-data/POSTPROC/${OLD_ARCHIVE} ${ONLINE_VALIDATION_DIR}/biofloat_ms
cp /marconi_scratch/usera07ogs/a07ogs01/Archive_biofloats_ms_validation_V2.tar   ${ONLINE_VALIDATION_DIR}/biofloat_ms
tar -xf ${ONLINE_VALIDATION_DIR}/biofloat_ms/${OLD_ARCHIVE} -C  ${ONLINE_VALIDATION_DIR}/biofloat_ms/V2
PREVIOUS_ARCH=${ONLINE_VALIDATION_DIR}/biofloat_ms/V2

#opa_prex_or_die "python biofloats_ms_plotter.py  -o ${ONLINE_VALIDATION_DIR}/biofloat_ms -p $PREVIOUS_ARCH -v ${ONLINE_VALIDATION_DIR} -a $OPA_ARCDIR_ROOT -d $OPA_RUNDATE "
#python biofloats_ms_plotter.py  -o ${ONLINE_VALIDATION_DIR_INPUT}/biofloat_ms -p $PREVIOUS_ARCH -v ${ONLINE_VALIDATION_DIR} -a $ARCHIVE_DIR -d $STARTTIME_s
echo "python biofloats_ms_plotter.py  -o ${ONLINE_VALIDATION_DIR}/biofloat_ms -p $PREVIOUS_ARCH -v ${ONLINE_VALIDATION_DIR_INPUT}/ANALYSIS_PROFILES_1week/ -a $ARCHIVE_DIR -d $OPA_RUNDATE"  #$STARTTIME_s"
python biofloats_ms_plotter.py  -o ${ONLINE_VALIDATION_DIR}/biofloat_ms -p $PREVIOUS_ARCH -v ${ONLINE_VALIDATION_DIR_INPUT}/ANALYSIS_PROFILES_1week/ -a $ARCHIVE_DIR -d $OPA_RUNDATE #$STARTTIME_s


exit 0

if [ 1 == 0 ]; then
mkdir -p $ONLINE_VALIDATION_DIR/PREVIOUS
mkdir -p $ONLINE_VALIDATION_DIR/ACTUAL

INPDIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/inpdir
SAT_WEEKLY_DIR=${INPDIR}/ONLINE/SAT/MULTISENSOR/1Km/NRT/WEEKLY_2_24/
SAT_DAILY_DIR=${INPDIR}/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/CHECKED_24/
export ONLINE_REPO=${INPDIR}/ONLINE


python archive_extractor.py --type analysis -st ${STARTTIME_a} -et ${END__TIME_a}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_s} -et ${STARTTIME_s}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS

fi
######


OPA_WRKDIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/wrkdir/2 #???
PREVIOUS_TUE_RUNDIR=${STARTTIME_s}
TMP_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/output/
TMP_DIR__ACT=$OPA_WRKDIR/MODEL/AVE_FREQ_1/

if [ 1 == 0 ]; then
f0_name=${ONLINE_VALIDATION_DIR}/Validation_f0_${PREVIOUS_TUE_RUNDIR}_on_weekly_Sat.${STARTTIME_s}.nc
opa_prex_or_die "python SatValidation_24.py -d ${STARTTIME_s} -f $TMP_DIR_PREV -s ${SAT_WEEKLY_DIR} -o $f0_name -m $MASKFILE" # MISFIT

for fc_day in 0 1 2; do
   DAY=$(date -d "$STARTTIME_s + $fc_day days"  +%Y%m%d )
   f_name=${ONLINE_VALIDATION_DIR}/Validation_f${fc_day}_${PREVIOUS_TUE_RUNDIR}_on_daily_Sat.${DAY}.nc
   opa_prex_or_die "python SatValidation_24.py  -d $DAY -f $TMP_DIR_PREV -s ${SAT_DAILY_DIR} -o ${f_name}  -m $MASKFILE"     # ERROR
done

# Questi devono essere archiviati in $ARCHIVE_DIR/POSTPROC/AVE_FREQ_1/validation/Sat ???

a0_name=${ONLINE_VALIDATION_DIR}/Validation_a0_${OPA_RUNDATE}_on_weekly_Sat.${STARTTIME_s}.nc
opa_prex_or_die "python SatValidation_24.py -d ${STARTTIME_s} -f $TMP_DIR__ACT -s ${SAT_WEEKLY_DIR} -o $a0_name  -m $MASKFILE"  # MISFIT

for ac_day in 1 2; do
   DAY=$(date -d "$STARTTIME_s + $ac_day days"  +%Y%m%d )
   a_name=${ONLINE_VALIDATION_DIR}/Validation_a${ac_day}_${OPA_RUNDATE}_on_daily_Sat.${DAY}.nc
   opa_prex_or_die "python SatValidation_24.py  -d $DAY -f $TMP_DIR__ACT -s ${SAT_DAILY_DIR} -o ${a_name}  -m $MASKFILE"    # ERROR
done

######
fi


PROFILERDIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/PROFILATORE
IMG_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/IMG
    BIO_DIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/output

PROFILERDIRA=${ONLINE_VALIDATION_DIR}/ACTUAL/PROFILATORE
IMG_DIR__ACT=${ONLINE_VALIDATION_DIR}/ACTUAL/IMG
    BIO_DIRA=${WRKDIR}/MODEL/AVE_FREQ_1


# Matchup analysis
mkdir -p $IMG_DIR_PREV
mkdir -p $IMG_DIR__ACT
#${STARTTIME_a} -et ${END__TIME_a} confronto Float-Analysis previous
#${STARTTIME_f} -et ${END__TIME_f} confronto Float con {analisi attuale, IMG_DIR__ACT}, {forecast, IMG_DIR_PREV}
python float_extractor.py -st ${STARTTIME_a} -et ${END__TIME_a} -i ${BIO_DIRP}  -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
#Matchup forecasts
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i ${BIO_DIRP}  -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
# Matchup actual
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i ${BIO_DIRA}  -b $PROFILERDIRA  -o $IMG_DIR__ACT -m $MASKFILE

mkdir -p ${ONLINE_VALIDATION_DIR}/matchup_outputs
python profileplotter_3.py -p $IMG_DIR_PREV -a $IMG_DIR__ACT -o  ${ONLINE_VALIDATION_DIR}/matchup_outputs  -f  ${ONLINE_VALIDATION_DIR}/BioFloats_Descriptor.xml # scorre tutti i files e genera le immagini a 3

######### biofloats matchups validation #######

mkdir -p ${ONLINE_VALIDATION_DIR}/ANALYSIS_over_2_weeks/TMP  ${ONLINE_VALIDATION_DIR}/ANALYSIS_over_2_weeks/PROFILATORE/PROFILES
for I in `ls ${TMP_DIR_PREV}/*nc ` ; do
  ln -fs $I ${ONLINE_VALIDATION_DIR}/ANALYSIS_over_2_weeks/TMP
done
for I in `ls ${TMP_DIR__ACT}/*nc ` ; do
  ln -fs $I ${ONLINE_VALIDATION_DIR}/ANALYSIS_over_2_weeks/TMP
done

for I in `ls ${ONLINE_VALIDATION_DIR}/PREVIOUS/PROFILATORE/PROFILES/*nc ` ; do
  ln -fs $I ${ONLINE_VALIDATION_DIR}/ANALYSIS_over_2_weeks/PROFILATORE/PROFILES
done
for I in `ls ${ONLINE_VALIDATION_DIR}/ACTUAL/PROFILATORE/PROFILES/*nc  ` ; do
  ln -fs $I ${ONLINE_VALIDATION_DIR}/ANALYSIS_over_2_weeks/PROFILATORE/PROFILES
done


OLD_ARCHIVE=Archive_biofloats_ms_validation_V4.tar
mkdir -p ${ONLINE_VALIDATION_DIR}/biofloat_ms/V4
cp ${OPA_ETCDIR}/static-data/POSTPROC/${OLD_ARCHIVE} ${ONLINE_VALIDATION_DIR}/biofloat_ms
tar -xf ${ONLINE_VALIDATION_DIR}/biofloat_ms/${OLD_ARCHIVE} -C  ${ONLINE_VALIDATION_DIR}/biofloat_ms/V4
PREVIOUS_ARCH=${ONLINE_VALIDATION_DIR}/biofloat_ms/V4

opa_prex_or_die "python biofloats_ms_plotter.py  -o ${ONLINE_VALIDATION_DIR}/biofloat_ms -p $PREVIOUS_ARCH -v ${ONLINE_VALIDATION_DIR} -a $OPA_ARCDIR_ROOT -d $OPA_RUNDATE "

###################################
exit 0


# PARALLEL MAP GENERATION -- done in phase_C1


BACKGROUND=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/etc/static-data/POSTPROC/background_medeaf.png
MAPS=/marconi_scratch/usera07ogs/a07ogs01/MAPS/
mkdir -p $MAPS/ORIG
cp /marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/etc/static-data/POSTPROC/fonts/TitilliumWeb-Regular.ttf $MAPS
cd $MAPS
BITSEA=/marconi/home/usera07ogs/a07ogs01/MAPPE/bit.sea/

# native
MODELDIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/wrkdir/2/MODEL/AVE_FREQ_1
XML_FILE=/marconi/home/usera07ogs/a07ogs01/MAPPE/bit.sea/postproc/Plotlist_bio.xml
mpirun -np 36 python $BITSEA/build_layer_maps.py -b $BACKGROUND -o $MAPS/ORIG -m $MASKFILE -i $MODELDIR -p $XML_FILE -g ave*N1p.nc
echo NATIVE DONE

# phys
MODELDIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/wrkdir/2/MODEL/AVE_PHYS
XML_FILE=/marconi/home/usera07ogs/a07ogs01/MAPPE/bit.sea/postproc/Plotlist_phys.xml
mpirun -np 36 python $BITSEA/build_layer_maps.py -b $BACKGROUND -o $MAPS/ORIG -m $MASKFILE -i $MODELDIR -p $XML_FILE -g ave*phys.nc
echo PHYS DONE


rm -rf $MAPS/OUT
mkdir $MAPS/OUT
cd $BITSEA/validation/online
BINDIR=/marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/HOST/marconi/bin/
mpirun -np 36 python map_compressor.py -i $MAPS/ORIG -o $MAPS/OUT -b $BINDIR



