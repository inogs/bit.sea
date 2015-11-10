import os, glob
from instruments.index_reader import index_reader
import instruments.mooring 

INPUTDIR      = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/mooring/"
OUTPUTDIR      = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/mooring_ave/"


INDEXED=index_reader()

fileLIST = glob.glob(INPUTDIR + "*.nc")

for filename in fileLIST:
    fileout = OUTPUTDIR + os.path.basename(filename)
    M = instruments.mooring.Mooring.fromfile(filename)
    varlist = M.available_params.rstrip(' ')
    for ivar, var in enumerate(varlist):
        if var == '': varlist.pop(ivar)
    for var in varlist:
        PRES, VAR = M.read(var,None)
        idepth = 0
        
    
    


