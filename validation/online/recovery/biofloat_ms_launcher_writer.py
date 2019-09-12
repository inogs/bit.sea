from commons.Timelist import TimeList, TimeInterval
from datetime import datetime
from dateutil.relativedelta import relativedelta

TI = TimeInterval("20190701","20190915","%Y%m%d")
DIR="/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/AVE/ANALYSIS"

TL = TimeList.fromfilenames(TI, DIR, "*nc", filtervar="N3n")

WEEKLY_REQS=TL.getWeeklyList(2)
maskfile="/gpfs/work/OGS_prod_0/OPA/V5C/prod/wrkdir/2/MODEL/meshmask.nc"
PROFILATORE="/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/PROFILATORE"

print "#! /bin/bash"

dirname = "/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/FLOAT_WEEKLY_VALIDATION_ARCHIVE/"
print "mkdir -p " + dirname

for req in WEEKLY_REQS:
    d = datetime.strptime(req.string, "%Y%m%d")
    #next_week = d+relativedelta(days=7)
    
    NRT_biofloat_out_name = "%s/BioFloat_Weekly_validation_%s.nc"  %(dirname,req.string )
    print "python biofloats_ms.py -d %s -m %s -b %s -o %s" %(req.string, maskfile, PROFILATORE, NRT_biofloat_out_name)