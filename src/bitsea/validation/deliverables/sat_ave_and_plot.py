import argparse
from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates
    Map_SATchl_{year}.png
    Maplog_SATchl_{year}.png
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = existing_dir_path,
                            required =True,
                            help = ''' Output dir'''
                            )

    parser.add_argument(   '--inputdir', '-i',
                            type = existing_dir_path,
                            required = True,
                            help = 'Input dir')
    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
                                required = True,
                                help = 'Path of the mask file')
    parser.add_argument(   '--starttime','-s',
                                type = date_from_str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-e',
                                type = date_from_str,
                                required = True,
                                help = 'start date in yyyymmdd format')


    return parser.parse_args()
args = argument()

import numpy as np
import matplotlib
matplotlib.use('Agg')
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask

from bitsea.layer_integral.mapplot import mapplot
from bitsea.layer_integral.mapplot import mapplotlog
import bitsea.commons.timerequestors as requestors
import bitsea.Sat.SatManager as Sat
import matplotlib.pyplot as pl
from bitsea.layer_integral import coastline

clon,clat = coastline.get()
TheMask = Mask.from_file(args.maskfile)
_,jpj,jpi = TheMask.shape

INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir



TI = TimeInterval.fromdatetimes(args.starttime,args.endtime)
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc", prefix="", dateformat="%Y%m%d")
if args.endtime.year > args.starttime.year:
    req_label=f"{args.starttime.year}_{args.endtime.year}"
else:
    req_label=f"{args.starttime.year}"


req = requestors.Generic_req(TI)
indexes,weights = TL.select(req)
nFrames = len(indexes)
SAT_3D=np.zeros((nFrames,jpj,jpi), np.float32)

for iFrame, k in enumerate(indexes):
    t = TL.Timelist[k]
    inputfile = TL.filelist[k]
    CHL = Sat.readfromfile(inputfile,'CHL')
    SAT_3D[iFrame,:,:] = CHL

Sat2d=Sat.averager(SAT_3D)

mask=TheMask.mask_at_level(0)
Sat2d[Sat2d<0] = np.nan
Sat2d[~mask] = np.nan
var = 'SATchl'

outfile    = OUTPUTDIR / f"Map_{var}_{req_label}.png"
print(outfile)
fig,ax     = mapplot({'clim':[0,0.4],  'data':Sat2d},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat, colormap='viridis')
ax.set_xlim([-6,36])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
#ax.ticklabel_format(fontsize=10)
ax.tick_params(axis='x', labelsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
ax.xaxis.set_ticks(np.arange(-6,36,4))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % (req_label, var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)

outfile    = OUTPUTDIR / f"Maplog_{var}_{req_label}.png"
print(outfile)
fig,ax     = mapplotlog({'clim':[1.01,1],  'data':Sat2d},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat, colormap='viridis')

ax.set_xlim([-6,36])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
#ax.ticklabel_format(fontsize=10)
ax.tick_params(axis='x', labelsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
ax.xaxis.set_ticks(np.arange(-6,16,6))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % (req_label , var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)
