#! /bin/bash

#export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V7C/
export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V7C/
echo $ONLINE_REPO

DATE_start=20181117
DATE_end=20211118
VARNAME='O2o'   # O2o:DOX
DATA_DAY=20211117  #  DATE_end -1
OUTDIR=/g100_scratch/userexternal/camadio0/SUPERFLOAT

python superfloat_o2o.py -Tst $DATE_start -Tend $DATE_end -Tday $DATA_DAY -o $OUTDIR -v $VARNAME

python dump_index.py -i $OUTDIR -o ${OUTDIR}/Float_Index.txt -t superfloat
echo ${OUTDIR}/Float_Index.txt
