import numpy as np
import scipy.io.netcdf as NC

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from layer_integral.mapbuilder import MapBuilder
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.layer import Layer


def TimeAverager(Filelist,weights,varname,mask):
    n=len(Filelist)
    jpk, jpj, jpi= mask.shape
    MSUM=np.zeros((jpk,jpj,jpi),np.float32)
    for t in range(n):
        filename=Filelist[t]
        #De      = DataExtractor(TheMask,filename,varname)
        #M = De.values
        ncIN = NC.netcdf_file(filename,"r")
        M = ncIN.variables[varname].data[0,:,:,:].copy()
        ncIN.close()
        MSUM += M*weights[t]
    averaged = MSUM/weights.sum()
    return averaged


def NCwriter(M2d,varname,outfile,mask):
    ncOUT = NC.netcdf_file(outfile,'w')
    _, jpj, jpi= mask.shape
    ncOUT.createDimension("longitude", jpi)
    ncOUT.createDimension("latitude", jpj)
    ncvar = ncOUT.createVariable(varname, 'f', ('latitude','longitude'))
    ncvar[:] = M2d
    ncOUT.close()

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')

INPUTDIR  = "/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/"
OUTPUTDIR = "/pico/scratch/usera07ogs/a07ogs02/GB/"
# Limito la richiesta, e poi prendo tutto quello che viene
#TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc", 'postproc/IOnames.xml')

VARLIST=['DIC','AC_','pH_','pCO']



TI = TimeInterval('20140404','20150629',"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc", 'postproc/IOnames.xml')
layer = Layer(0,200)
Seas_reqs = TL.getSeasonList()

for req in Seas_reqs:
    print req
    indexes, weights = TL.select(req)
    prefix = req.string.replace(" ","")
    for var in VARLIST[:1]:
        filelist=[]
        for k in indexes:
            t = TL.Timelist[k]
            filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
            filelist.append(filename)
        M3d = TimeAverager(filelist, weights, var, TheMask)
        De      = DataExtractor(TheMask,rawdata=M3d)
        integrated = MapBuilder.get_layer_average(De, layer)
        outfile    = OUTPUTDIR + prefix + '.' +  var + ".nc"




# Now, the year
import commons.timerequestors as requestors
MY_YEAR = TimeInterval('20140401','20150401',"%Y%m%d")
req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)

for var in VARLIST[:1]:
    filelist=[]
    for k in indexes:
        t = TL.Timelist[k]
        filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
        filelist.append(filename)
    M3d     = TimeAverager(filelist, weights, var, TheMask)
    De      = DataExtractor(TheMask,rawdata=M3d)
    integrated = MapBuilder.get_layer_average(De, layer)
    outfile    = OUTPUTDIR + "year." + var + "nc"

