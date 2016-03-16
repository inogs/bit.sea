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
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output dir'''
                            )

    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required = True,
                            default = './',
                            help = 'Input dir')


    return parser.parse_args()



def NCwriter(M2d,varname,outfile,mask):
    ncOUT = NC.netcdf_file(outfile,'w')
    _, jpj, jpi= mask.shape
    ncOUT.createDimension("longitude", jpi)
    ncOUT.createDimension("latitude", jpj)
    ncvar = ncOUT.createVariable(varname, 'f', ('latitude','longitude'))
    ncvar[:] = M2d
    ncOUT.close()

args = argument()
coast=np.load('Coastline.npy')
clon=coast['Lon']
clat=coast['Lat']
TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')

# ANALYSIS AND FORECAST PRE OPERATIONAL QUALIFICATION RUN
#INPUTDIR  = "/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/"
#OUTPUTDIR = "/pico/home/userexternal/gcossari/COPERNICUS/Carbonatic17/MAPPE_MEDIE/"

# REANALYSIS V2 RUN
INPUTDIR  = args.inputdir#"/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/MODEL/AVE_FREQ_2/"
OUTPUTDIR = args.outdir#"/pico/home/userexternal/gcossari/COPERNICUS/REANALYSIS_V2/MAPPE_MEDIE/"

LIMIT_PER_MASK=[5,5,5]

#VARLIST=['DIC','AC_','PH_','pCO']
LAYERLIST=[Layer(0,10),Layer(0,50),Layer(0,150)]
VARLIST=['ppn','N1p','N3n','PH_','pCO','P_l'] # saved as mg/m3/d --> * Heigh * 365/1000
VARUNI=['gC/m^2/y','mmol/m^3','mmol/m^3','','ppm','mmol/m^3'];
CLIM=[[0, 200],[0, 0.1],[0, 4],[7.9, 8.2],[300,480],[0, 1]];
VARCONV=[365./1000.,1,1,1,1,1]
MEDIA_O_INTEGRALE=[1,0,0,0,0,0]
#MEDIA_O_INTEGRALE=1 # 1 -> INTEGRALE:  * heigth of the layer
                    # 0 -> MEDIA    :  average of layer
TI = TimeInterval('20000101','20121230',"%Y%m%d") # VALID FOR REANALYSIS RUN
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*N1p.nc", 'postproc/IOnames.xml')

# CHOICE OF THE TIME SELECTION
import commons.timerequestors as requestors

MY_YEAR = TimeInterval('20000101','20121230',"%Y%m%d") # requestor generico per la media del reanalysis 1999-2012
req_label='Ave:1999-2014'

req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)

for iv, var in enumerate(VARLIST[:1]):
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

        if MEDIA_O_INTEGRALE[iv]==1:
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
            integrated=integrated * Hlayer * VARCONV[iv]
        else:
            integrated=integrated * VARCONV[iv]

#        mask200=TheMask.mask_at_level(200)
        mask200=TheMask.mask_at_level(LIMIT_PER_MASK[il])
#        clim = [M3d[TheMask.mask].min(), M3d[TheMask.mask].max()]
        clim=CLIM[iv]
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

        #pl.show(block=False)
        pl.close(fig)
        #NCwriter(integrated,var,outfile,TheMask)
