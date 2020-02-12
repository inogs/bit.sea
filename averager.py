import numpy as np

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer
from commons import netcdf3
from layer_integral.mapbuilder import MapBuilder
from layer_integral.mapplot import *
from commons.dataextractor import DataExtractor
from commons.time_averagers import TimeAverager3D
import matplotlib.pyplot as pl
from commons import season


TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')

INPUTDIR  = "/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/"
OUTPUTDIR = "/pico/scratch/userexternal/gcossari/Carbonatic-17/"

LAYERLIST=[Layer(0,50), Layer(50,100), Layer(100,150), Layer(150,200), Layer(200,500), Layer(500,1000), Layer(1000,1500), Layer(1500,4000)]
VARLIST=['DIC','AC_','PH_','pCO']



TI = TimeInterval('20140404','20150629',"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc")
s = season.season()
s.setseasons(["1221","0321","0622","0921"])
Seas_reqs = TL.getSeasonList(s)


# AC_2014.aut.0000-0050m.nc
for req in Seas_reqs:
    print req
    indexes, weights = TL.select(req)
    prefix = req.string.replace(" ","")
    for var in VARLIST:
        print var
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
            clim = [M3d[TheMask.mask].min(), M3d[TheMask.mask].max()]
            fig,ax     = mapplot({'varname':var, 'clim':clim, 'layer':layer, 'data':integrated, 'date':req.string},fig=None,ax=None,mask=TheMask)
            outfile    = OUTPUTDIR + prefix + '.' +  var + "." + layer.longname() +  ".nc"
            netcdf3.write_2d_file(integrated,var,outfile,TheMask)
            outfile    = OUTPUTDIR + var + "." + prefix  +  "." + layer.longname() +  ".png"
            fig.savefig(outfile)
            pl.close(fig)







# Now, the whole year
import commons.timerequestors as requestors
MY_YEAR = TimeInterval('20140401','20150401',"%Y%m%d")
req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)

for var in VARLIST:
    print var
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
        clim = [M3d[TheMask.mask].min(), M3d[TheMask.mask].max()]
        fig,ax     = mapplot({'varname':var, 'clim':clim, 'layer':layer, 'data':integrated, 'date':'annual'},fig=None,ax=None,mask=TheMask)
        outfile    = OUTPUTDIR + var + "." + "annual."  + layer.longname() + ".png"
        fig.savefig(outfile)
        pl.close(fig)
