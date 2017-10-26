#! /bin/bash

#PBS -N MAPS
#PBS -l walltime=0:30:00
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=100GB
#PBS -q route
#PBS -A OGS_dev_0

# cd $PBS_O_WORKDIR
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


OPA_RUNDATE=20171017 # tuesday 

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


mkdir -p $ONLINE_VALIDATION_DIR/PREVIOUS
mkdir -p $ONLINE_VALIDATION_DIR/ACTUAL

if [ 1 == 0 ]; then

python archive_extractor.py --type analysis -st ${STARTTIME_a} -et ${END__TIME_a}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_s} -et ${STARTTIME_s}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS
python archive_extractor.py --type forecast -st ${STARTTIME_f} -et ${END__TIME_f}  -a ${ARCHIVE_DIR}  -o ${ONLINE_VALIDATION_DIR}/PREVIOUS


mkdir -p $ONLINE_VALIDATION_DIR/PREVIOUS/output_phys
mkdir -p $ONLINE_VALIDATION_DIR/ACTUAL/output_phys

#Cutting the physics domain INGV to OGS:
python ingv_cutter.py -i ${ONLINE_VALIDATION_DIR}/PREVIOUS/output_phys_ingv -o ${ONLINE_VALIDATION_DIR}/PREVIOUS/output_phys -M $INGV_MASK -m $MASKFILE
python ingv_cutter.py -i /marconi/home/usera07ogs/a07ogs00/OPA/V3C-dev/wrkdir/2/OPAOPER   -o ${ONLINE_VALIDATION_DIR}/ACTUAL/output_phys -M $INGV_MASK -m $MASKFILE -l *T.nc -p "" -f "%y%m%d"
if [ $? -ne 0 ] ; then echo "ERROR" ; exit 1 ; fi



PROFILERDIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/PROFILATORE
IMG_DIR_PREV=${ONLINE_VALIDATION_DIR}/PREVIOUS/IMG
   PHYS_DIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/output_phys
    BIO_DIRP=${ONLINE_VALIDATION_DIR}/PREVIOUS/output_bio

PROFILERDIRA=${ONLINE_VALIDATION_DIR}/ACTUAL/PROFILATORE
IMG_DIR__ACT=${ONLINE_VALIDATION_DIR}/ACTUAL/IMG
   PHYS_DIRA=${ONLINE_VALIDATION_DIR}/ACTUAL/output_phys
    BIO_DIRA=${WRKDIR}/MODEL/AVE_FREQ_1

# Matchup analysis
mkdir -p $IMG_DIR_PREV
mkdir -p $IMG_DIR__ACT
#${STARTTIME_a} -et ${END__TIME_a} confronto Float-Analysis previous
#${STARTTIME_f} -et ${END__TIME_f} confronto Float con {analisi attuale, IMG_DIR__ACT}, {forecast, IMG_DIR_PREV}
python float_extractor.py -st ${STARTTIME_a} -et ${END__TIME_a} -i ${BIO_DIRP} -p ${PHYS_DIRP}  -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
#Matchup forecasts
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i ${BIO_DIRP} -p ${PHYS_DIRP}  -b $PROFILERDIRP  -o $IMG_DIR_PREV -m $MASKFILE
# Matchup actual
python float_extractor.py -st ${STARTTIME_f} -et ${END__TIME_f} -i ${BIO_DIRA} -p ${PHYS_DIRA}  -b $PROFILERDIRA  -o $IMG_DIR__ACT -m $MASKFILE

mkdir -p ${ONLINE_VALIDATION_DIR}/matchup_outputs
python profileplotter_3.py -p $IMG_DIR_PREV -a $IMG_DIR__ACT -o  ${ONLINE_VALIDATION_DIR}/matchup_outputs  -f  ${ONLINE_VALIDATION_DIR}/BioFloats_Descriptor.xml # scorre tutti i files e genera le immagini a 3
fi


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
MODELDIR=/marconi_scratch/usera07ogs/a07ogs01/online_validation_data/ACTUAL/output_phys
XML_FILE=/marconi/home/usera07ogs/a07ogs01/MAPPE/bit.sea/postproc/Plotlist_phys.xml
mpirun -np 36 python $BITSEA/build_layer_maps.py -b $BACKGROUND -o $MAPS/ORIG -m $MASKFILE -i $MODELDIR -p $XML_FILE -g ave*votemper.nc
echo PHYS DONE

rm -rf $MAPS/OUT
mkdir $MAPS/OUT
cd $BITSEA/validation/online
mpirun -np 36 python map_compressor.py -i $MAPS/ORIG -o $MAPS/OUT



