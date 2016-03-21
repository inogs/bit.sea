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
import Sat.SatManager as Sat
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


INPUTDIR  = args.inputdir#"/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/MODEL/AVE_FREQ_2/"
OUTPUTDIR = args.outdir#"/pico/home/userexternal/gcossari/COPERNICUS/REANALYSIS_V2/MAPPE_MEDIE/"



TI = TimeInterval('20000101','20121230',"%Y%m%d") # VALID FOR REANALYSIS RUN
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc", 'postproc/IOnames.xml')

# CHOICE OF THE TIME SELECTION
import commons.timerequestors as requestors

MY_YEAR = TimeInterval('20000101','20121230',"%Y%m%d") # requestor generico per la media del reanalysis 1999-2012
req_label='Ave:1999-2014'

req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)
nFrames = len(indexes)
SAT_3D=np.zeros((nFrames,jpj,jpi), np.float32)

### READ SATELLITE
filelist=[]
for k in indexes:
    t = TL.Timelist[k]
    filename = INPUTDIRSAT + t.strftime("%Y%m") + "_d-OC_CNR-L4-CHL-MedOC4_SAM_7KM-MED-REP-v02.nc"
    CHL = Sat.readfromfile(inputfile)
    SAT_3D[iFrame,:,:] = CHL

Sat2d=Sat.averager(SAT_3D)

Sat2d = TimeAverager2DSAT(filelist, weights, "lchlm", TheMask)
Sat2d[Sat2d>1.e19]=np.nan
Sat2d=Sat2d*mask200
Sat2d[Sat2d==0]=np.nan


    fig,ax     = mapplot({'varname':var, 'clim':clim, 'layer':layer, 'data':integrated200, 'date':'annual'},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
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
