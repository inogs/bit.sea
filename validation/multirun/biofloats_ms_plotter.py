import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates
     - in the first output directory time series png files similar to them of mdeeaf web site.
     - in the 2nd output directory   two tables,  _BIAS.txt and _RMSE.txt for each variable
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = '')

    parser.add_argument(   '--figdir', '-f',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--tabledir', '-t',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()


from commons.time_interval import TimeInterval
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import numpy as np
from commons.layer import Layer
from basins import V2 as OGS
from commons.utils import addsep
# from profiler import TL
import scipy.io.netcdf as NC
from commons.utils import writetable
from datetime import datetime
from profiler_2015 import *

OUT_FIGDIR        = addsep(args.figdir)
OUT_TABLEDIR       = addsep(args.tabledir)
inputfile        = args.inputfile

class ncreader():
    def __init__(self, filename):
        ncIN = NC.netcdf_file(filename,"r")
        self.nVAR = ncIN.dimensions['var']
        self.nSUB = ncIN.dimensions['sub']
        self.nDEPTH  = ncIN.dimensions['depth']
        self.SUBLIST = ncIN.sublist.split(",")
        self.LAYERLIST=ncIN.layerlist.split(",")
        self.VARLIST = ncIN.varlist.split(",")
        self.npoints = ncIN.variables['npoints'].data.copy()
        self.bias    = ncIN.variables['bias'   ].data.copy()
        self.rmse    = ncIN.variables['rmse'   ].data.copy()
    def plotdata(self,VAR,var,sub,depth):
        ivar = self.VARLIST.index(var)
        isub = self.SUBLIST.index(sub)
        idepth= self.LAYERLIST.index(depth)
        return VAR[ivar,:, isub,idepth]



DATAfile = ncreader(inputfile)

def single_plot(longvar, var, sub, layer, timeinterval ):
    if TL.inputFrequency == 'weekly':
        times = TL.Timelist
    else:
        if TL.inputFrequency=='daily':
             WEEKS=TL.getWeeklyList(5)
             times = [datetime(w.year,w.month, w.day)  for w in WEEKS]
    bias1 = DATAfile.plotdata(DATAfile.bias   , var, sub, layer)
    rmse1 = DATAfile.plotdata(DATAfile.rmse   , var, sub, layer)
    numb1 = DATAfile.plotdata(DATAfile.npoints, var, sub, layer)
    
    ii=numb1==0
    rmse1[ii] = np.nan
    bias1[ii] = np.nan

    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()

    ax2.bar(times,numb1,width=7, color='0.5', alpha=0.3, align='center')
    ax2.set_ylabel(' n. of BGC-Argo floats', fontsize=20)
    ax2.set_ylim([0,numb1.max() +2])
#    ax2.set_yticklabels(ax2.yaxis.get_major_ticks(),fontsize=16)


    ax.plot(times,bias1,'m.-', label='bias')
    ax.plot(times,rmse1,'k.-', label='rmse')
    if longvar == 'Chlorophyll' : 
	ax.set_ylim([-0.4, 0.4])
	ax.set_ylabel('bias, rmse mg/m$^3$', fontsize=20)
#	ax.ticklabel_format(axis='both', fontsize=20)

    if longvar == 'Nitrate'     : 
	ax.set_ylim([-5, 5])
	ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=20)
    if longvar == 'Oxygen'      :
	ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=20)
    #if longvar == 'Oxygen'      : ax.set_ylim([-40, 40])
        
#    ax.set_ylabel('bias, rmse mg/m$^3$')
    ax.legend(loc=2)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for l in ax2.get_yticklabels() : l.set_fontsize(16)

    ii = np.zeros((len(times),) , np.bool)
    for k,t in enumerate(times) : ii[k] = timeinterval.contains(t)
    biasm = np.nanmean(bias1[ii])
    rmsem = np.nanmean(rmse1[ii])
    return fig, biasm, rmsem, ax, ax2

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
VARLIST = ['Chla','N3n','O2o']
VARLONGNAMES=['Chlorophyll','Nitrate','Oxygen']
VARLIST = ['P_l','N3n']
VARLONGNAMES=['Chlorophyll','Nitrate']
SUBLIST = OGS.NRT3.basin_list
nSub = len(SUBLIST)
nLayers = len(LAYERLIST)
#ti_restrict = TimeInterval("20150101","20170101","%Y%m%d")
ti_restrict = TimeInterval(DATESTART,DATE__END,"%Y%m%d")

column_names=[layer.string() for layer in LAYERLIST]
row_names   =[sub.name for sub in SUBLIST]

for ivar, var in enumerate(VARLIST):
    BIAS = np.zeros((nSub,nLayers),np.float32)
    RMSE = np.zeros((nSub,nLayers),np.float32)
    for isub, sub in enumerate(SUBLIST):
        for ilayer, layer in enumerate(LAYERLIST):
            outfile = "%s%s.%s.%s.png" % (OUT_FIGDIR,var,sub.name,layer.longname())
            print outfile
            fig,bias,rmse,ax,ax2  = single_plot(VARLONGNAMES[ivar],var,sub.name,layer.string(), ti_restrict)
            BIAS[isub,ilayer] = bias
            RMSE[isub,ilayer] = rmse
            title = "%s %s %s " %(VARLONGNAMES[ivar], sub.extended_name, layer.string())
            fig.suptitle(title, fontsize=20)
            fig.savefig(outfile)
            pl.close(fig)

    writetable(OUT_TABLEDIR +  var + '_BIAS.txt',BIAS,row_names, column_names)
    writetable(OUT_TABLEDIR +  var + '_RMSE.txt',RMSE,row_names, column_names)
    
