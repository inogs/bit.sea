# limstd and dayf are information on how the maps used to calculate sat variance have been obtained:
## dayf is the number of days over which satellite mapsa are averaged
## limstd is the number of std used to exclude outlyers with respect to climatology
# These values are not used in var_sat but to assign name of output files

INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/WEEKLY_24/
OUTVAR_SAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/VAR_SAT/
DAYF=7 # see file-beginning note
LIMSTD=2 # see file-beginning note
MESH=Mesh24
#MESH=V4mesh

mkdir -p $OUTVAR_SAT


echo mpirun -np 12 python var_sat.py -i $INSAT -o $OUTVAR_SAT -m $MESH -s $LIMSTD -f $DAYF
#More than 12 procs are unuseful becasue of iteration on 12 climatological months


INMOD=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/CHLA_HC24/
OUTVAR_DIFF=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/VAR_DIFF_MOD_SAT/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc

mkdir -p $OUTVAR_DIFF

echo mpirun -np 12 python VarianceDiffModObs.py -s $INSAT -d $INMOD -o $OUTVAR_DIFF -m $MESH -f $MASKFILE
#More than 12 procs are unuseful becasue of iteration on 12 climatological months


OUTVAR_MOD=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/VAR_MOD/

mkdir -p $OUTVAR_MOD

echo python varmod_vdiffvsat.py -s $OUTVAR_SAT -d $OUTVAR_DIFF -o $OUTVAR_MOD -m $MASKFILE -v 0.5 0.75
