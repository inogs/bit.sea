#! /bin/bash
for filename in `ls /gpfs/scratch/userexternal/eterzic0/KD_MONTHLY/V4/*nc `; do 
   I=`basename $filename`
   localname=ave.${I:0:6}15-00:00:00.KD490.nc
   ln -fs $filename $localname
done
