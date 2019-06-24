DIRCHL=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/CFR_CLIMA/DA_chl/
mkdir -p DA_chllinkall

for dir in $(ls -d $DIRCHL*); do
   echo $dir
   for filechl in $(ls $dir/chl.*); do
       echo $filechl
       ln -s $filechl DA_chllinkall/
   done
done
