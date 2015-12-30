import numpy as np
import scipy.io.netcdf as NC

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer

from layer_integral.mapbuilder import MapBuilder
from commons.dataextractor import DataExtractor
from commons.time_averagers import TimeAverager3D

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
OUTPUTDIR = "/pico/scratch/userexternal/gbolzon0/GB/"

LAYERLIST=[Layer(0,50), Layer(50,100), Layer(100,150), Layer(150,200), Layer(200,500), Layer(500,1000), Layer(1000,1500), Layer(1500,4000)]
VARLIST=['DIC','AC_','pH_','pCO']



TI = TimeInterval('20140404','20150629',"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc", 'postproc/IOnames.xml')
Seas_reqs = TL.getSeasonList()

for req in Seas_reqs:
    print req
    indexes, weights = TL.select(req)
    prefix = req.string.replace(" ","")
    for var in VARLIST:
        # setting up filelist for requested season -----------------
        filelist=[]
        for k in indexes:
            t = TL.Timelist[k]
            filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
            filelist.append(filename)
        # ----------------------------------------------------------

        M3d = TimeAverager3D(filelist, weights, var, TheMask)
        for layer in LAYERLIST:
            print layer
            De         = DataExtractor(TheMask,rawdata=M3d)
            integrated = MapBuilder.get_layer_average(De, layer)
            outfile    = OUTPUTDIR + prefix + '.' +  var + ".nc"
            NCwriter(integrated,var,outfile,TheMask)







# Now, the whole year
import commons.timerequestors as requestors
MY_YEAR = TimeInterval('20140401','20150401',"%Y%m%d")
req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)

for var in VARLIST:
    # setting up filelist for requested period -----------------
    filelist=[]
    for k in indexes:
        t = TL.Timelist[k]
        filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
        filelist.append(filename)
    # ----------------------------------------------------------
    M3d     = TimeAverager3D(filelist, weights, var, TheMask)
    for layer in LAYERLIST:
        De      = DataExtractor(TheMask,rawdata=M3d)
        integrated = MapBuilder.get_layer_average(De, layer)
        outfile    = OUTPUTDIR + "year." + var + "nc"
        NCwriter(integrated,var,outfile,TheMask)

