
export MASKFILE=/pico/scratch/userexternal/lmariott/FDA_all2015_newstd_3days/wrkdir/MODEL/meshmask.nc

#LOC=/pico/scratch/userexternal/gbolzon0/MASSIMILI/FLOAT_DA_00/RESULTS
LOC=/pico/scratch/userexternal/gbolzon0/MASSIMILI/FDA_all2015_newstd_3days/RESULTS
mkdir -p $LOC/Dens_N1p $LOC/Dens_N3n $LOC/Dens_chl

mkdir -p $LOC/validation_tables/
#python profiler.py

python massimili_valid.py  -o $LOC/validation_tables/ -m $MASKFILE

python density_plots.py     -M $MASKFILE -o $LOC/Dens_N1p  -v N1p -m 0
python density_plots.py     -M $MASKFILE -o $LOC/Dens_N3n  -v N3n -m 0
python density_plots.py     -M $MASKFILE -o $LOC/Dens_chl  -v P_l -m 0



