import argparse
# python sat_model_RMSD_and_plot.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE  -o Fig4.1/ -s 20170101 -e 20180101

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

    parser.add_argument(   '--satdir', '-s',
                            type = str,
                            required = True,
                            default = './',
                            help = 'Satellite dir')

    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required = True,
                            default = './',
                            help = 'Input dir')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')
 
    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-et',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')


    return parser.parse_args()
args = argument()

import numpy as np
import os
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer
from commons.dataextractor import DataExtractor

from layer_integral.mapplot import mapplot
from layer_integral.mapplot import mapplotlog
import commons.timerequestors as requestors
import Sat.SatManager as Sat
import matplotlib.pyplot as pl
from layer_integral import coastline
from layer_integral.mapbuilder import MapBuilder
import matchup.matchup as matchup

clon,clat = coastline.get()
TheMask=Mask(args.maskfile)
_,jpj,jpi = TheMask.shape

MODEL_DIR  = args.inputdir
SATDIR    = args.satdir
OUTPUTDIR = args.outdir

START_TIME=args.starttime
END___TIME=args.endtime

surf_layer = Layer(0,10)
mask0_2D = TheMask.mask_at_level(0.0)
coastmask = mask0_2D

dateformat ="%Y%m%d"
#TI = TimeInterval('20170101','20180101',"%Y%m%d") 
TI = TimeInterval(START_TIME,END___TIME,"%Y%m%d")
sat_TL = TimeList.fromfilenames(TI, SATDIR,"*.nc", prefix="", dateformat="%Y%m%d")
print sat_TL
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*P_l.nc")

suffix = os.path.basename(sat_TL.filelist[0])[8:]

#MY_YEAR = TimeInterval('20170101','20180101',"%Y%m%d")
MY_YEAR = TimeInterval(START_TIME,END___TIME,"%Y%m%d")
req_label='RMSD:2019'

req = requestors.Generic_req(MY_YEAR)
indexes,weights = model_TL.select(req)
nFrames = len(indexes)
SD=np.zeros((nFrames,jpj,jpi), np.float32)
COUNTS=np.zeros((jpj,jpi), np.float32)
tmp_COUNTS=np.zeros((nFrames,jpj,jpi), np.float32)

for iFrame, k in enumerate(indexes):
    modeltime = model_TL.Timelist[k]

    CoupledList = sat_TL.couple_with([modeltime])
    sattime = CoupledList[0][0]
    satfile = SATDIR + "/" + sattime.strftime(dateformat) + suffix
    modfile = model_TL.filelist[k]
    print modfile

    De         = DataExtractor(TheMask,filename=modfile, varname='P_l')
    Model      = MapBuilder.get_layer_average(De, surf_layer)
    print modeltime
    Sat24 = Sat.readfromfile(satfile,var='CHL') 

    cloudsLand = (np.isnan(Sat24)) | (Sat24 > 1.e19) | (Sat24<0)
    modelLand  = np.isnan(Model) #lands are nan
    nodata     = cloudsLand | modelLand
    selection = ~nodata & coastmask
    Model[~selection] = np.nan
    Sat24[~selection] = np.nan
    tmp_COUNTS[iFrame,:,:][selection]=1
    M= matchup.matchup(Model,Sat24)
 
    SD[iFrame,:,:][selection] = ((M.Model-M.Ref)**2)
    SD[iFrame,:,:][~selection]=np.nan

COUNTS=np.sum(tmp_COUNTS,axis=0)
Sat2d=np.sqrt(np.nanmean(SD,axis=0))

mask=TheMask.mask_at_level(0)
Sat2d[Sat2d<0] = np.nan
Sat2d[~mask] = np.NaN
var = 'SATchl'

fig,ax     = mapplot({'clim':[0.001,0.05],  'data':Sat2d},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
ax.set_xlim([-5,36])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
ax.ticklabel_format(fontsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
#ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label +  ".png"
ax.xaxis.set_ticks(np.arange(-2,36,6))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % ('annual RMSD', var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)


fig,ax     = mapplotlog({'clim':[0.01,1],  'data':Sat2d},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
ax.set_xlim([-5,36])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
ax.ticklabel_format(fontsize=10)
ax.text(-4,44.5,var + ' [mg /m^3]',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
#ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
outfile    = OUTPUTDIR + "Maplog_" + var + "_" + req_label +  ".png"
ax.xaxis.set_ticks(np.arange(-2,36,6))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % ('annual', var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)

fig,ax     = mapplot({'clim':[40,55],  'data':COUNTS},fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
ax.set_xlim([-5,36])
ax.set_ylim([30,46])
ax.set_xlabel('Lon').set_fontsize(12)
ax.set_ylabel('Lat').set_fontsize(12)
ax.ticklabel_format(fontsize=10)
ax.text(-4,44.5,var + ': # POINTS ',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
#ax.text(-4,32,'Ave:' + layer.string() ,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
outfile    = OUTPUTDIR + "Map_COUNTS_" + var + "_" + req_label +  ".png"
ax.xaxis.set_ticks(np.arange(-2,36,6))
ax.yaxis.set_ticks(np.arange(30,46,4))
ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
ax.grid()
title = "%s %s" % ('Number of records ', var)
fig.suptitle(title)
fig.savefig(outfile)
pl.close(fig)
