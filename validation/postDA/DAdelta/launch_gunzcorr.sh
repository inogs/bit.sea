DIRRST=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/DA_RST/

for dir in $(ls -d $DIRRST*); do
   echo $dir
   for filecorr in $(ls $dir/*_corr*gz); do
       echo $filecorr
       gunzip -f $filecorr
   done
done
