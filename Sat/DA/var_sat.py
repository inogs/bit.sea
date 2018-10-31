from commons import timerequestors
from commons.Timelist import TimeInterval, TimeList
from Sat import SatManager
from postproc import masks
import numpy as np
from layer_integral.mapplot import mapplot
import sys

maskSat = getattr(masks,"Mesh24")
jpi = maskSat.jpi
jpj = maskSat.jpj
INPUTDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/CCI1km_Interp24_10days/"
OUTPUTDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/VARSAT/"

TL=TimeList.fromfilenames(None, INPUTDIR, "*.nc", prefix="", dateformat="%Y%m%d")
valCax=[.4,.4,.9,1,.6,.2,.1,.1,.1,.3,.5,.6]

for month in range(1,13):
    req = timerequestors.Clim_month(month)
    ii,w = TL.select(req)
    filelist=[TL.filelist[k] for k in ii if TL.Timelist[k].month==month] # we consider strictly the month

    chl2DT =np.zeros((jpj, jpi), np.float32)
    chl2D2T=np.zeros((jpj, jpi), np.float32)
    NPoints=np.zeros((jpj, jpi), np.int32)
    
    for filename in filelist:
        A = SatManager.readfromfile(filename)
        bad = np.isnan(A) | (A==-999)
        A[bad] = 0
        app = (~bad).astype(np.int32)
        
        chl2DT += A
        chl2D2T += A**2
        NPoints += app

    noPoints=NPoints==0
    NPoints[noPoints] = 1 # just to divide
    chl2Dm  = chl2DT/NPoints
    chl2D2m = chl2D2T/NPoints
    var2D   = chl2D2m - chl2Dm**2
    varFig=np.sqrt(var2D)
    varFig[noPoints] = np.nan    
    var2D[noPoints] = -999.0

    NPoints[noPoints] = -999 # just to save on file

    outfile_var= "%svar2Dsat.CCI.F10.2.%02d.nc"  %( OUTPUTDIR, month)
    outfile__np=       "%snp.CCI.F10.2.%02d.nc"  %( OUTPUTDIR, month)
    outfile_png=   "%sdevstd.CCI.F10.2.%02d.png" %( OUTPUTDIR, month)
    SatManager.dumpGenericNativefile(outfile_var, var2D, "variance", maskSat)
    SatManager.dumpGenericNativefile(outfile__np, NPoints, "np"    , maskSat)

    
    fig, ax = mapplot({'data':varFig, 'clim':[0,valCax[month-1]]})
    ax.set_title(req.longname())
    fig.savefig(outfile_png)

