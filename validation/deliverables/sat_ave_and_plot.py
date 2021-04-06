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
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')


    return parser.parse_args()
args = argument()

import numpy as np

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
#from commons.layer import Layer

from layer_integral.mapplot import mapplot
from layer_integral.mapplot import mapplotlog
import commons.timerequestors as requestors
import Sat.SatManager as Sat
import matplotlib.pyplot as pl
from layer_integral import coastline

clon,clat = coastline.get()
TheMask=Mask(args.maskfile)
_,jpj,jpi = TheMask.shape

INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir



TI = TimeInterval('20190101','20200101',"%Y%m%d") # VALID FOR REANALYSIS RUN
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc", prefix="", dateformat="%Y%m%d")
print TL.Timelist


#MY_YEAR = TimeInterval('20190101','20200101',"%Y%m%d") 
req_label='2019_2020'
#req_label='201901'

#req = requestors.Generic_req(MY_YEAR)
req = requestors.Generic_req(TI)
indexes,weights = TL.select(req)
nFrames = len(indexes)
SAT_3D=np.zeros((nFrames,jpj,jpi), np.float32)

for iFrame, k in enumerate(indexes):
    t = TL.Timelist[k]
    inputfile = TL.filelist[k]
    CHL = Sat.readfromfile(inputfile,'CHL')
#    CHL = Sat.readfromfile(inputfile,'lchlm')
    SAT_3D[iFrame,:,:] = CHL

Sat2d=Sat.averager(SAT_3D)

mask=TheMask.mask_at_level(0)
Sat2d[Sat2d<0] = np.nan
Sat2d[~mask] = np.NaN
var = 'SATchl'
#layer = Layer(0,10)

#fig,ax     = mapplot({'varname':var, 'clim':[0,0.4], 'layer':layer, 'data':Sat2d, 'date':'annual'},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
fig,ax     = mapplot({'clim':[0,0.4],  'data':Sat2d},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
ax.set_xlim([-9,16])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
ax.ticklabel_format(fontsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
#ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label +  ".png"
ax.xaxis.set_ticks(np.arange(-9,16,4))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % ('annual', var)
title = "%s %s" % (req_label, var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)


fig,ax     = mapplotlog({'clim':[0.01,1],  'data':Sat2d},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
ax.set_xlim([-9,16])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
ax.ticklabel_format(fontsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
#ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
outfile    = OUTPUTDIR + "Maplog_" + var + "_" + req_label +  ".png"
ax.xaxis.set_ticks(np.arange(-9,16,6))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % ('annual', var)
#title = "%s %s" % (req_label , var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)
