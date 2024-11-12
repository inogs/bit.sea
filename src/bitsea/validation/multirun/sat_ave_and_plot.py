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
import os

from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer

from bitsea.layer_integral.mapplot import mapplot
import bitsea.commons.timerequestors as requestors
import bitsea.Sat.SatManager as Sat
import matplotlib.pyplot as pl
from bitsea.layer_integral import coastline

clon,clat = coastline.get()
maskfile = os.getenv('MASKFILE')
TheMask = Mask.from_file(maskfile)
_,jpj,jpi = TheMask.shape

Timestart=os.getenv("START_DATE")
Time__end=os.getenv("END_DATE")


INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir



TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc", prefix="", dateformat="%Y%m%d")


#MY_YEAR = TimeInterval('20130101','20140101',"%Y%m%d") 
req_label = Timestart[0:4] #'Ave:2013'

req = requestors.Generic_req(TI)
indexes,weights = TL.select(req)
nFrames = len(indexes)
SAT_3D=np.zeros((nFrames,jpj,jpi), np.float32)

for iFrame, k in enumerate(indexes):
    t = TL.Timelist[k]
    inputfile = INPUTDIR + t.strftime("%Y%m%d") + "_d-OC_CNR-L4-CHL-MedOC4_SAM_7KM-MED-REP-v02.nc"
    CHL = Sat.readfromfile(inputfile,'lchlm')
    SAT_3D[iFrame,:,:] = CHL

Sat2d=Sat.averager(SAT_3D)

masknan=TheMask.mask_at_level(0)
Sat2d[~masknan] = np.NaN
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
