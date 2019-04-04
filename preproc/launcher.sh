#! /bin/bash

INPUTDIR=/gpfs/work/OGS_prod_0/OPA/V5C/devel/wrkdir/2/MODEL/FORCINGS/ 
OUTPUTDIR=DeltaT/
INGVMASK=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc


mpirun -np 17 python delta_t_from_uvw.py -i $INPUTDIR -o $OUTPUTDIR -m $INGVMASK

python decision_deltat.py -i $OUTPUTDIR

