import numpy as np
from commons.Timelist import TimeInterval,TimeList
import os
import Sat.SatManager as Sat
from postproc.masks import V4mesh

Bias_monthly_file="/galileo/home/userexternal/ateruzzi/SATBIAS_CCI_MULTI/biasMonthPoint.npy"
INPUTDIR="/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST2019/WEEKLY_16/"
# INPUTDIR="/marconi_scratch/userexternal/ateruzzi/MULTIsat_biascorrect/WEEKLY_16/"
# INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/DT/WEEKLY/"
# OUTPUTDIR="/pico/scratch/userexternal/gbolzon0/REANALYSIS_2016/wrkdir/SAT_CORRECTED/"
# CHECKDIR="/pico/scratch/userexternal/gbolzon0/REANALYSIS_2016/wrkdir/SAT_CHECK/"
# OUTPUTDIR="/marconi_scratch/userexternal/ateruzzi/MULTIsat_biascorrect/SAT_CORRECTED/"
OUTPUTDIR="/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST2019/SAT_CORRECTED/"
CHECKDIR="/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST2019/SAT_CHECK/"
threshold = 0.01

TI = TimeInterval("20190101","20200101","%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR, "*nc", prefix="", dateformat="%Y%m%d")


BIAS=np.load(Bias_monthly_file)
ii = np.isnan(BIAS)
BIAS[ii] = -999.0
for iFrame in range(12):
    outfile = CHECKDIR + "bias%02d.nc" %(iFrame+1)
    bias=BIAS[:,:,iFrame]
    Sat.dumpGenericNativefile(outfile, bias, 'lchlm', V4mesh)



BIAS[ii]=0.0 # null correction on no-data points

for iFrame, filename in enumerate(TL.filelist):
    outfile=OUTPUTDIR + os.path.basename(filename)
    Sat_orig = Sat.readfromfile(filename,var="CHL")
    goods = Sat_orig > 0
    index_month = TL.Timelist[iFrame].month -1
    print outfile, index_month
    bias = BIAS[:,:,index_month]
    Sat_corrected = np.ones_like(Sat_orig, dtype=np.float32)*(-999.0)
    Sat_corrected[goods] = Sat_orig[goods] + bias[goods]
    
    negative = (Sat_orig > 0) &  (Sat_corrected<= 0)
    print "negative points : ", negative.sum()
    Sat_corrected[negative] = threshold
    Sat.dumpGenericNativefile(outfile, Sat_corrected, "lchlm", V4mesh)
    

