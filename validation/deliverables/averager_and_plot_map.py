#
# ATTENZIONE  this script needs: py_env and pyload


import numpy as np
import scipy.io.netcdf as NC

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer

from layer_integral.mapbuilder import MapBuilder
from layer_integral.mapplot import *
from commons.dataextractor import DataExtractor
from commons.time_averagers import TimeAverager3D
import pylab as pl

def NCwriter(M2d,varname,outfile,mask):
    ncOUT = NC.netcdf_file(outfile,'w')
    _, jpj, jpi= mask.shape
    ncOUT.createDimension("longitude", jpi)
    ncOUT.createDimension("latitude", jpj)
    ncvar = ncOUT.createVariable(varname, 'f', ('latitude','longitude'))
    ncvar[:] = M2d
    ncOUT.close()

coast=np.load('Coastline.npy')
clon=coast['Lon']
clat=coast['Lat']
TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')

# ANALYSIS AND FORECAST PRE OPERATIONAL QUALIFICATION RUN
INPUTDIR  = "/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/"
OUTPUTDIR = "/pico/home/userexternal/gcossari/COPERNICUS/Carbonatic17/MAPPE_MEDIE/"

# REANALYSIS V2 RUN
INPUTDIR  = "/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA/wrkdir/MODEL/AVE_FREQ_2/"
OUTPUTDIR = "/pico/home/userexternal/gcossari/COPERNICUS/REANALYSIS_V2/MAPPE_MEDIE/"



#LAYERLIST=[Layer(0,50), Layer(50,100), Layer(100,150), Layer(150,200), Layer(200,500), Layer(500,1000), Layer(1000,1500), Layer(1500,4000)]
LAYERLIST=[Layer(0,200)]
#LAYERLIST=[Layer(40,60)]
#LAYERLIST=[Layer(0,50), Layer(180,220)]
#LAYERLIST=[Layer(0,50),Layer(0,150)]
#LIMIT_PER_MASK=[50,150]

#LAYERLIST=[Layer(0,10)]
LIMIT_PER_MASK=[5]

VARLIST=['DIC','AC_','PH_','pCO']

VARLIST=['ppn'] # saved as mg/m3/d --> * Heigh * 365/1000
VARUNI=['gC/m^2/y']; CLIM=[0, 200]; VARCONV=365./1000.
#VARLIST=['N1p'];VARUNI=['mmol/m^3'];CLIM=[0, 0.1]; VARCONV=1.
#VARLIST=['N3n'];VARUNI=['mmol/m^3'];CLIM=[0, 4]; VARCONV=1.
#VARLIST=['PH_'];VARUNI=[''];CLIM=[7.9, 8.2]; VARCONV=1.
#VARLIST=['pCO'];VARUNI=['ppm'];CLIM=[300,480]; VARCONV=1.
#VARLIST=['P_l'];VARUNI=['mg/m^3'];CLIM=[0, 1]; VARCONV=1.

MEDIA_O_INTEGRALE=1 # 1 -> INTEGRALE:  * heigth of the layer
                    # 0 -> MEDIA    :  average of layer  
#TI = TimeInterval('20140404','20150629',"%Y%m%d") # VALID FOR PRE-OPERATIONAL QUALIFICATION RUN
TI = TimeInterval('20000101','20121230',"%Y%m%d") # VALID FOR REANALYSIS RUN
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc", 'postproc/IOnames.xml')

# CHOICE OF THE TIME SELECTION
import commons.timerequestors as requestors
#MY_YEAR = TimeInterval('20140630','20150629',"%Y%m%d") # requestor generico per la media annuale da 7/14 a 6/15
#req_label='Ave:07/14-06/15'
MY_YEAR = TimeInterval('20000101','20121230',"%Y%m%d") # requestor generico per la media del reanalysis 1999-2012
req_label='Ave:1999-2014'

req = requestors.Generic_req(MY_YEAR)
#req = requestors.Clim_month(2); req_label='clim_Feb' # requestor for climatologia of a specific month
#req = requestors.Clim_month(5); req_label='clim_May' # requestor for climatologia of a specific month
#req = requestors.Clim_month(8); req_label='clim_Aug' # requestor for climatologia of a specific month
#req = requestors.Clim_month(11); req_label='clim_Nov' # requestor for climatologia of a specific month
indexes,weights = TL.select(req)

for iv, var in enumerate(VARLIST):
    print var
    varuni=VARUNI[iv]
    # setting up filelist for requested period -----------------
    filelist=[]
    for k in indexes:
        t = TL.Timelist[k]
        filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
        filelist.append(filename)
    # ----------------------------------------------------------
    M3d     = TimeAverager3D(filelist, weights, var, TheMask)
    for il,layer in enumerate(LAYERLIST):
        De      = DataExtractor(TheMask,rawdata=M3d)
        integrated = MapBuilder.get_layer_average(De, layer)

        if MEDIA_O_INTEGRALE==1:
# calcolo l'altezza del layer
            top_index = np.where(De._mask.zlevels >= layer.top)[0][0]
            bottom_index = np.where(De._mask.zlevels < layer.bottom)[0][-1]
        #Workaround for Python ranges
            bottom_index += 1
        #Build local mask matrix
            lmask = np.array(De._mask.mask[top_index:bottom_index,:,:], dtype=np.double)
        #Build dz matrix
            dzm = np.ones_like(lmask, dtype=np.double)
            j = 0
            for i in range(top_index, bottom_index):
                dzm[j,:,:] = De._mask.dz[i]
                j += 1
        #Build height matrix (2D)
            Hlayer = (dzm * lmask).sum(axis=0)
            integrated=integrated * Hlayer * VARCONV
        else:
            integrated=integrated * VARCONV

#        mask200=TheMask.mask_at_level(200)
        mask200=TheMask.mask_at_level(LIMIT_PER_MASK[il])
#        clim = [M3d[TheMask.mask].min(), M3d[TheMask.mask].max()]
        clim=CLIM
        integrated200=integrated*mask200 # taglio il costiero
        integrated200[integrated200==0]=np.nan # sostituisco gli 0 con i NAN

#change the colormap 
        pl.set_cmap('gray_r')
        fig,ax     = mapplot({'varname':var, 'clim':clim, 'layer':layer, 'data':integrated200, 'date':'annual'},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
        ax.set_xlim([-5,36])
        ax.set_ylim([30,46])
        ax.set_xlabel('Lon').set_fontsize(12)
        ax.set_ylabel('Lat').set_fontsize(12)
        ax.ticklabel_format(fontsize=10)
        ax.text(-4,44.5,var + ' [' + varuni + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
        if  MEDIA_O_INTEGRALE==1:
          ax.text(-4,32,'Int:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
          outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Int" + layer.longname() + ".png"
          outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Int_GRAY" + layer.longname() + ".png"
        else:
          ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
          outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Ave" + layer.longname() + ".png"
          outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Ave_GRAY" + layer.longname() + ".png"
        ax.xaxis.set_ticks(np.arange(-2,36,6))
        ax.yaxis.set_ticks(np.arange(30,46,4))
        ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')

        outfile=outfile.replace(":","").replace("/14","_14").replace("/15","_15")

        fig.savefig(outfile)

        pl.show(block=False)
        #pl.close(fig)
        #NCwriter(integrated,var,outfile,TheMask)

