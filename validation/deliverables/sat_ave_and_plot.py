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

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer

from layer_integral import mapplot
import commons.timerequestors as requestors
import Sat.SatManager as Sat
import pylab as pl

coast=np.load('Coastline.npy')
clon=coast['Lon']
clat=coast['Lat']
TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')
_,jpj,jpi = TheMask.shape()

INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir



TI = TimeInterval('20000101','20121230',"%Y%m%d") # VALID FOR REANALYSIS RUN
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc", 'postproc/IOnames.xml')


MY_YEAR = TimeInterval('20000101','20121230',"%Y%m%d") 
req_label='Ave:1999-2014'

req = requestors.Generic_req(MY_YEAR)
indexes,weights = TL.select(req)
nFrames = len(indexes)
SAT_3D=np.zeros((nFrames,jpj,jpi), np.float32)

for iFrame, k in enumerate(indexes):
    t = TL.Timelist[k]
    inputfile = INPUTDIR + t.strftime("%Y%m") + "_d-OC_CNR-L4-CHL-MedOC4_SAM_7KM-MED-REP-v02.nc"
    CHL = Sat.readfromfile(inputfile)
    SAT_3D[iFrame,:,:] = CHL

Sat2d=Sat.averager(SAT_3D)

var = 'SATchl'
layer = Layer(0,10)

fig,ax     = mapplot({'varname':var, 'clim':[0,0.4], 'layer':layer, 'data':Sat2d, 'date':'annual'},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
ax.set_xlim([-5,36])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
ax.ticklabel_format(fontsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Ave" + layer.longname() + ".png"
ax.xaxis.set_ticks(np.arange(-2,36,6))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
fig.savefig(outfile)
pl.close(fig)
