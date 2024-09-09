from commons import genUserDateList as DL
from commons.Timelist import TimeList

OUTDIR="/marconi_scratch/userexternal/gbolzon0/GLOBAL_REANALYSIS_BIO_001_018/ORIG/"
ACCESS="--user MED_OGS_TRIESTE_IT --pwd NEdifupa --motu http://my.cmems-du.eu/motu-web/Motu"
product="GLOBAL_REANALYSIS_BIO_001_018-TDS"

BINDIR="/marconi_scratch/userexternal/gbolzon0/motu-client-python-1.5.00-20180223190259664/src/python"
BOX=" --longitude-min -15 --longitude-max 0 --latitude-min 30 --latitude-max 47  --depth-min 0. --depth-max 5000. "


#PRODUCT=" --service-id GLOBAL_REANALYSIS_BIO_001_018-TDS  --product-id dataset-global-nahindcast-bio-001-018-V5-chl --variable CHL --variable nav_lon --variable nav_lat "
#TIME=" --date-min \"2016-12-01 12:00:00\" --date-max \"2016-12-31 12:00:00\"  "
#IO=" --out-dir " + OUTDIR +  "--out-name GLO_chl.20161201.nc "


months=DL.getTimeList("20150101-00:00:00", "20150501-00:00:00", "months=1")
TL=TimeList(months)

DATASET_DICT={'CHL': "dataset-global-nahindcast-bio-001-018-V5-chl",
              'OXY': "dataset-global-nahindcast-bio-001-018-V5-O2",}
VARS = DATASET_DICT.keys()

LINES=[]
LINES.append("#! /bin/bash \n")
LINES.append("BINDIR="+BINDIR + "\n")

MONTHLY_REQS=TL.getMonthlist()
for var in VARS:
    for req in MONTHLY_REQS:
        outfile = "ave.%s01-00:00:00.%s.nc"  %(req.string, var)
        dataset=DATASET_DICT[var]
        startime=req.time_interval.start_time.strftime("%Y-%m-%d %H:%M:%S")
        end__ime=req.time_interval.end_time.strftime("%Y-%m-%d %H:%M:%S")
        PRODUCT=" --service-id %s  --product-id %s --variable %s  --variable nav_lon --variable nav_lat " %(product,dataset,var)
        TIME=" --date-min \"%s\" --date-max \"%s\"  " %(startime, end__ime)
        IO = " --out-dir %s --out-name %s " %(OUTDIR, outfile)
        
        line = "python $BINDIR/motu-client.py  "  + ACCESS + PRODUCT + BOX + TIME + IO
        LINES.append(line + "\n")

F=open("motu_launcher.sh",'w')
F.writelines(LINES)
F.close()
