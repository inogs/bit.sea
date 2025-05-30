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

    parser.add_argument(   '--optype', '-t',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['integral','mean'],
                                help ="  INTEGRALE:  * heigth of the layer, MEDIA    :  average of layer")

    parser.add_argument(   '--layertype', '-l',
                                type = str,
                                required = True,
                                choices = ['multilayer','0-200layer'],
                                help ="  multilayer:  layers for validation, 0-200layer    :  layer 0-200 only")

    return parser.parse_args()
args = argument()

import numpy as np
import os

from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer

from bitsea.layer_integral.mapbuilder import MapBuilder
from bitsea.layer_integral.mapplot import mapplot,pl
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.time_averagers import TimeAverager3D
from bitsea.layer_integral import coastline
import bitsea.commons.timerequestors as requestors


clon,clat = coastline.get()

maskfile = os.getenv('MASKFILE')
TheMask = Mask.from_file(maskfile)

Timestart=os.getenv("START_DATE")
Time__end=os.getenv("END_DATE")



INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir
var       = args.varname

if args.layertype=='multilayer':
   LAYERLIST=[Layer(  0,  10), \
              Layer( 10,  50), \
              Layer( 50, 100), \
              Layer(100, 150), \
              Layer(150, 300), \
              Layer(300, 600), \
              Layer(600,1000)]
else:
    LAYERLIST=[Layer(0,200)]
#VARLIST=['ppn','N1p','N3n','PH_','pCO','P_l'] # saved as mg/m3/d --> * Heigh * 365/1000 #VARLIST=['DIC','AC_','PH_','pCO']
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

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc")


#MY_YEAR = TimeInterval('20130101','20140101',"%Y%m%d") 
req_label = Timestart[0:4] #'Ave:2013'

req = requestors.Generic_req(TI)
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

    if args.optype=='integral':
        integrated = MapBuilder.get_layer_integral(De, layer)
    else:
        integrated = MapBuilder.get_layer_average(De, layer)  
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
