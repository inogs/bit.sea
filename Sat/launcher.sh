#! /bin/bash

PYPATH=/pico/home/usera07ogs/a07ogs00/OPA/V2C/HOST/pico/lib/python2.7/site-packages/
export PYTHONPATH=$PYTHONPATH:$PYPATH
 
#WEEKLY_1KMDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_1km_Friday/
#WEEKLY_16DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_16_Friday/
#MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc
#mpirun -np 10 python interpolator.py -i $WEEKLY_1KMDIR -o $WEEKLY_16DIR -inmesh SAT1km_mesh -m V4mesh -M $MASKFILE


CHECKED_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/CHECKED/
WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_1km
#mpirun python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR -m SAT1km_mesh -t weekly_tuesday



#multisensor 1km
CHECKED_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/CHECKED/
WEEKLY_DIR=/pico/scratch/userexternal/ateruzzi/SATELLITE/CCI16_1km/WEEKLY_1km/
#mpirun -np 10 python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR -m SAT1km_mesh -t weekly_tuesday


WEEKLY_1KMDIR=$WEEKLY_DIR
WEEKLY_16_DIR=/pico/scratch/userexternal/ateruzzi/SATELLITE/CCI16_1km/WEEKLY_16/
MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc

# Multisensor DT
WEEKLY_1KMDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/DT/WEEKLY_2_1km/
WEEKLY_1KMDIR=/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST2018/WEEKLY_2_1km_2018
WEEKLY_16_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST2018/WEEKLY_2_16_2018/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc


mkdir -p $WEEKLY_16_DIR
echo mpirun -np 10 python interpolator.py -i $WEEKLY_1KMDIR -o $WEEKLY_16_DIR -inmesh SAT1km_mesh -m V4mesh -M $MASKFILE


exit 0


# Section SAT CHECK

ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG/
CHECKDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/CHECKED/
CLIM_FILE=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SatClimatology.nc

python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIM_FILE -m SAT1km_mesh

ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/ORIG/
CHECKDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/CHECKED/
python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIM_FILE -m SAT1km_mesh
