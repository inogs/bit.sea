import argparse

from commons.time_interval import TimeInterval
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import numpy as np
from commons.layer import Layer
from basins import V2 as OGS
from commons.utils import addsep
from profiler import TL
import scipy.io.netcdf as NC

OUTFIG_DIR       = "biofloats/" #addsep("args.outdir")

inputfile="pippo.nc"

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

def single_plot(longvar, var, sub, layer ):
    times = TL.Timelist
    bias1 = DATAfile.plotdata(DATAfile.bias   , var, sub, layer)
    rmse1 = DATAfile.plotdata(DATAfile.rmse   , var, sub, layer)
    numb1 = DATAfile.plotdata(DATAfile.npoints, var, sub, layer)
    
    ii=numb1==0
    rmse1[ii] = np.nan
    bias1[ii] = np.nan

    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()

    ax2.bar(times,numb1,width=7, color='0.5', alpha=0.3, align='center')
    ax2.set_ylabel(' # Points')
    ax2.set_ylim([0,numb1.max() +2])


    ax.plot(times,bias1,'m.-', label='bias')
    ax.plot(times,rmse1,'k.-', label='rmse')
    if longvar == 'Chlorophyll' : ax.set_ylim([-0.5, 0.5])
    if longvar == 'Nitrate'     : ax.set_ylim([-10, 10])
    if longvar == 'Oxygen'      : ax.set_ylim([-40, 40])
        
    ax.set_ylabel('bias, rmse mg/m$^3$')
    ax.legend(loc=2)
    ax2.legend(loc=1)
    ax.set_title(longvar)

    return fig

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
VARLONGNAMES=['Chlorophyll','Nitrate','Oxygen']
SUBLIST = OGS.NRT3.basin_list


for ivar, var in enumerate(VARLIST):
    for isub, sub in enumerate(SUBLIST[:2]):
        for layer in LAYERLIST[:2]:
            outfile = "%s%s.%s.%s.png" % (OUTFIG_DIR,var,sub.name,layer.longname())
            print outfile
            fig = single_plot(VARLONGNAMES[ivar],var,sub.name,layer.string())
            fig.savefig(outfile)
            pl.close(fig)
