import os, glob
from instruments.mooring import *
import datetime
from commons.Timelist import TimeList

INPUTDIR      = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/mooring/"
OUTPUTDIR     = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/mooring_ave/"




fileLIST = glob.glob(INPUTDIR + "*.nc")

for filename in fileLIST:
    fileout = OUTPUTDIR + os.path.basename(filename)
    M = Mooring.fromfile(filename)
    if M is None :
        print filename
        continue
    varlist = M.available_params.rsplit(' ')
    for ivar, var in enumerate(varlist):
        if var == '': varlist.pop(ivar)
        
    var = varlist[0]
    PRES, VAR, TIME = M.read_raw(var)
    datelist=[]
    for t in TIME:
        dateprofile = basetime + datetime.timedelta(days = t)
        datelist.append(dateprofile)
    TL = TimeList(datelist)
    import sys
    sys.exit()
#     for var in varlist:
#         PRES, VAR, TIME = M.read_raw(var)
#         idepth = 0
