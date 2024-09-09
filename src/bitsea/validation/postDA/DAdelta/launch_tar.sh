DIRTAR=/gss/gss_work/DRES_OGS_BiGe/ateruzzi/RA_COAST/2016/output/DA_DIRS/

for dir in $(ls $DIRTAR/2016*.tar); do
   filename=$(basename $dir)
   date=${filename:0:8}
   echo $date
   OUTDIR=DA_RST/$date
   mkdir -p $OUTDIR
   tar -xvf $dir RST*
   mv RST*$date*nc* $OUTDIR/
done


