DIRTAR=/gss/gss_work/DRES_OGS_BiGe/ateruzzi/RA_COAST/2017/output/DA_DIRS/

for dir in $(ls $DIRTAR/201*.tar); do
   filename=$(basename $dir)
   date=${filename:0:8}
   echo $date
   OUTDIR=DA_chl/$date
   mkdir -p $OUTDIR
   tar -xvf $dir chl*.nc*
   mv chl.${date}-12*nc* $OUTDIR/
done


