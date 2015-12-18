
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
INPUTDIR="/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/"
# Limito la richiesta, e poi prendo tutto quello che viene
TI = TimeInterval('20140401','20150401',"%Y%m%d")

TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc", 'postproc/IOnames.xml')

VARLIST=['DIC','AC_','pH_','pCO']

for var in VARLIST:
    for t in TL.Timelist:
        filelist=[]
        filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
        filelist.append(filename)
        
