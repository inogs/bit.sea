#! /bin/bash

PYPATH=/pico/home/usera07ogs/a07ogs00/OPA/V2C/HOST/pico/lib/python2.7/site-packages/
export PYTHONPATH=$PYTHONPATH:$PYPATH
 
CHECKED_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/CHECKED/
WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_1km
#mpirun python aveSat.py -i $CHECKED_DIR -o $WEEKLY_DIR -m SAT1km_mesh -t weekly_tuesday



#multisensor 1km
CHECKED_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/CHECKED/
MONTHLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/MONTHLY_1km
mpirun -np 10 python aveSat.py -i $CHECKED_DIR -o $MONTHLY_DIR -m SAT1km_mesh -t monthly


MONTHLY_1KMDIR=$MONTHLY_DIR
MONTHLY_24_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/MONTHLY_24
MASKFILE=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc
cd CHL
mpirun -np 10 python interpolator.py -i $MONTHLY_1KMDIR -o $MONTHLY_24_DIR -m Mesh24 -M $MASKFILE




# Section SAT CHECK

ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG/
CHECKDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/CHECKED/
CLIM_FILE=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SatClimatology.nc

python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIM_FILE -m SAT1km_mesh

ORIGDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/ORIG/
CHECKDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/CHECKED/
python sat_check.py -i $ORIGDIR -o $CHECKDIR -c $CLIM_FILE -m SAT1km_mesh
