#! /bin/bash

#export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V7C/
export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V7C/
echo $ONLINE_REPO

DATE_start=20181117
DATE_end=20211118
VARNAME='O2o'   # O2o:DOX
DATA_DAY=20211117  #  DATE_end -1

OUTDIR=/g100_scratch/userexternal/camadio0/SUPERFLOAT

if [ -d "$OUTDIR" ]; then
   echo "'$OUTDIR' found "
else
   mkdir $OUTDIR
fi

python superfloat_chla.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
python superfloat_o2o.py -Tst $DATE_start -Tend $DATE_end -Tday $DATE_DAY -o $OUTDIR -v $VARNAME
python superfloat_nitrate.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
python superfloat_par.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
python superfloat_bbp700.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
python superfloat_ph.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
python dump_index.py -i $OUTDIR -o ${OUTDIR}/Float_Index.txt -t superfloat

echo ${OUTDIR}/Float_Index.txt

Condition_to_copy=false
if [ "$Condition_to_copy" = true ] ; then
   cp LAUNCHER_SUPERFLOAT.sh $OUTDIR
fi

Concat_reports=false
if [ "$Concat_reports" = true ] ; then
   REP_DIR=/g100/home/userexternal/camadio0/float_preproc/TimeSeries_plots/Superfloat_oxy/OUTPUTS/
   python append_report.py -rep $REP_DIR
   echo $Concat_reports
fi







