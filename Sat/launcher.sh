#! /bin/bash


#multisensor 1km
CHECKED_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/DAILY/CHECKED/
MONTHLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/MONTHLY_1km
mpirun -np 10 python aveSat.py -i $CHECKED_DIR -o $MONTHLY_DIR -m SAT1km_mesh -t monthly