# limstd and dayf are information on how the maps used to calculate sat variance have been obtained:
## dayf is the number of days over which satellite mapsa are averaged
## limstd is the number of std used to exclude outlyers with respect to climatology
# These values are not used in var_sat but to assign name of output files

INSAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_SAT/WEEKLYweight/WEEKLY_24/
OUTVAR_SAT=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/VAR_SAT_CHECKED/
DAYF=7 # see file-beginning note
LIMSTD=2 # see file-beginning note
MESH=Mesh24
#MESH=V4mesh

mkdir -p $OUTVAR_SAT


echo mpirun -np 12 python var_sat.py -i $INSAT -o $OUTVAR_SAT -m $MESH -s $LIMSTD -f $DAYF
#More than 12 procs are unuseful becasue of iteration on 12 climatological months


#INMOD=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/CHLA_HC24/
#OUTVAR_DIFF=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/VAR_DIFF_MOD_SAT/
#MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc

#mkdir -p $OUTVAR_DIFF

#echo mpirun -np 12 python VarianceDiffModObs.py -s $INSAT -d $INMOD -o $OUTVAR_DIFF -m $MESH -f $MASKFILE
#More than 12 procs are unuseful becasue of iteration on 12 climatological months

INMOD2DAYS=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/CHLA_HC24/
MODWEEKLY=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/CHLA_HC24_WEEKLY/
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_REA/meshmask.nc

mkdir -p $MODWEEKLY

echo mpirun -np 32 python aveModWeighted.py -i $INMOD2DAYS -o $MODWEEKLY -m $MASKFILE -t weekly_monday -s 3.5



OUTVAR_MOD=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/PURE_VAR_MOD_WEEKLY/
OUTVAR_MOD=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/PURE_VAR_MOD_WEEKLYlimits/

mkdir -p $OUTVAR_MOD
echo mpirun -np 12 python varmod_ratiosat.py -s $OUTVAR_SAT -d $MODWEEKLY -o $OUTVAR_MOD -m $MASKFILE -v 0.5 0.75 -r 15 15

#OUTVAR_MOD=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/VAR_MOD/

#mkdir -p $OUTVAR_MOD

#echo python varmod_vdiffvsat.py -s $OUTVAR_SAT -d $OUTVAR_DIFF -o $OUTVAR_MOD -m $MASKFILE -v 0.5 0.75

INSTATP=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/postprocHC24/output/AVE_FREQ_1/STAT_PROFILES/
NPYEOF=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/EOF_PKL/

mkdir -p $NPYEOF

echo python eofs_statp.py -i $INSTATP -o $NPYEOF -m $MASKFILE -v Chla -d 200


NCEOF=/gpfs/scratch/userexternal/ateruzzi/REA_24_DA_static/EOF_NC/
NEOF=26
MAPCNPY=/galileo/home/userexternal/ateruzzi/INTERP_MAPSERkmeans/mapser.19771215.24.npy
GRID3DVAR=/gpfs/scratch/userexternal/ateruzzi/DA_staticTRANSITION_V6/DA_static_data/3D_VAR/GRID/BFM_grid__Sat.nc

mkdir -p $NCEOF


echo mpirun -np 12 python eofs_nc.py -i $NPYEOF -v $OUTVAR_MOD -n $NEOF -o $NCEOF -m $MASKFILE -c $MAPCNPY -g $GRID3DVAR




