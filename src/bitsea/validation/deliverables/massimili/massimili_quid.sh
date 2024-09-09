
export MASKFILE=/marconi_scratch/userexternal/ateruzzi/DA_FLOAT_SAT/RUN_FLOAT_2015_01/wrkdir/MODEL/meshmask.nc

LOC=OUTPUT
mkdir -p $LOC/Dens_N1p $LOC/Dens_N3n $LOC/Dens_chl $LOC/validation_tables/
#python profiler.py

python massimili_valid.py  -o $LOC/validation_tables/ -m $MASKFILE

python density_plots.py     -M $MASKFILE -o $LOC/Dens_N1p  -v N1p -m 0
python density_plots.py     -M $MASKFILE -o $LOC/Dens_N3n  -v N3n -m 0
python density_plots.py     -M $MASKFILE -o $LOC/Dens_chl  -v P_l -m 0



