DIRRST=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/DA_RST/

for dir in $(ls -d $DIRRST*); do
   echo $dir
   for filemis in $(ls $dir/*chl_mis*gz); do
       echo $filemis
       gunzip -f $filemis
   done
done
