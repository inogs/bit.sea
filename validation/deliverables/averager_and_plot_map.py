import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot something
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

    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['P_i','N1p', 'N3n', 'pCO','pH','ppn'] )

    parser.add_argument(   '--mapdepthfilter', '-m',
                                type = float,
                                required = False,
                                default = 0.0,
                                help  = 'The level chosen to filter coast data' )
    parser.add_argument(   '--top', 
                                type = float,
                                required = True,
                                help  = 'Top of the layer, for integration' )
    parser.add_argument(   '--bottom',
                                type = float,
                                required = True,
                                help  = 'Bottom of the layer, for integration' )

    parser.add_argument(   '--optype', '-t',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['integral','mean'],
                                help ="  INTEGRALE:  * heigth of the layer, MEDIA    :  average of layer")

    return parser.parse_args()
args = argument()

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


INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir
var       = args.varname
layer     = Layer(args.top, args.bottom)


LAYERLIST=[layer]
VARLIST=['ppn','N1p','N3n','PH_','pCO','P_l'] # saved as mg/m3/d --> * Heigh * 365/1000 #VARLIST=['DIC','AC_','PH_','pCO']
UNITS_DICT={
         'ppn' : 'gC/m^2/y',
         'N1p' : 'mmol /m^3',
         'N3n' : 'mmol /m^3',
         'PH'  : '',
         'pCO2': 'ppm',
         'P_i' :'mmol /m^3'
         }

CLIM_DICT={
         'ppn' : [0, 220],
         'N1p' : [0, 0.15],
         'N3n' : [0, 4],
         'PH'  : [7.9, 8.2],
         'pCO2': [300,480],
         'P_i' : [0, 0.4]
         }


CONVERSION_DICT={
         'ppn' : 365./1000,
         'N1p' : 1,
         'N3n' : 1,
         'PH'  : 1,
         'pCO2': 1,
         'P_i' : 1
         }

TI = TimeInterval('20000101','20121230',"%Y%m%d") # VALID FOR REANALYSIS RUN
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc", 'postproc/IOnames.xml')

# CHOICE OF THE TIME SELECTION
import commons.timerequestors as requestors

#MY_YEAR = TimeInterval('20000101','20121230',"%Y%m%d") # requestor generico per la media del reanalysis 1999-2012
MY_YEAR = TimeInterval('19990101','20141230',"%Y%m%d") 
req_label='Ave:1999-2014'

req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)

VARCONV=CONVERSION_DICT[var]
# setting up filelist for requested period -----------------
filelist=[]
for k in indexes:
    t = TL.Timelist[k]
    filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + ".nc"
    filelist.append(filename)
# ----------------------------------------------------------
M3d     = TimeAverager3D(filelist, weights, var, TheMask)
for il,layer in enumerate(LAYERLIST):
    De      = DataExtractor(TheMask,rawdata=M3d)
    integrated = MapBuilder.get_layer_average(De, layer)

    if args.optype=='integral':
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
    mask=TheMask.mask_at_level(args.mapdepthfilter)
#        clim = [M3d[TheMask.mask].min(), M3d[TheMask.mask].max()]
    clim=CLIM_DICT[var]
    integrated_masked=integrated*mask # taglio il costiero
    integrated_masked[integrated_masked==0]=np.nan # sostituisco gli 0 con i NAN


    #pl.set_cmap('gray_r') #changes the colormap
    fig,ax     = mapplot({'varname':var, 'clim':clim, 'layer':layer, 'data':integrated_masked, 'date':'annual'},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
    ax.set_xlim([-5,36])
    ax.set_ylim([30,46])
    ax.set_xlabel('Lon').set_fontsize(12)
    ax.set_ylabel('Lat').set_fontsize(12)
    ax.ticklabel_format(fontsize=10)
    ax.text(-4,44.5,var + ' [' + UNITS_DICT[var] + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
    if  args.optype=='integral':
        ax.text(-4,32,'Int:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
        outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Int" + layer.longname() + ".png"
    else:
        ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
        outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Ave" + layer.longname() + ".png"
    ax.xaxis.set_ticks(np.arange(-2,36,6))
    ax.yaxis.set_ticks(np.arange(30,46,4))
    ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
    ax.grid()
    fig.savefig(outfile)
    pl.close(fig)
