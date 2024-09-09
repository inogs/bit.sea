DIRCHL=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/CFR_CLIMA/DA_chl/

for dir in $(ls -d $DIRCHL*); do
   echo $dir
   for filechl in $(ls $dir/chl.*gz); do
       echo $filechl
       gunzip -f $filechl
   done
done
