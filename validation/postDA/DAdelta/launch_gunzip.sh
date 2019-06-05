DIRRST=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/DAdelta/DA_RST/

for dir in $(ls -d $DIRRST/201*); do
   echo $dir
   for filerst in $(ls $dir/RST*gz); do
       echo $filerst
       gunzip -f $filerst
   done
done
