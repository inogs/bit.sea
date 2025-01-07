import argparse
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates
     - time series png files similar to them of mdeeaf web site.
     - 4 tables,for each variable
       var_BIAS.txt
       var_RMSE.txt
       var_REF.txt
       var_MOD.txt
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputfile','-i',
                                type = existing_file_path,
                                required = True)

    parser.add_argument(   '--outdir', '-o',
                                type = existing_dir_path,
                                required = True)
    parser.add_argument(   '--basedir', '-b',
                                type = existing_dir_path,
                                required = True,
                                help = """ PROFILATORE dir, already generated""")
    parser.add_argument(   '--var', '-v',
                                type = str,
                                choices = ['P_l','N3n','O2o','P_c','POC'],
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()


from bitsea.commons.Timelist import TimeList, TimeInterval
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import numpy as np
from bitsea.commons.layer import Layer
from bitsea.basins import V2 as OGS
import netCDF4 as NC

from bitsea.commons.utils import writetable, nanmean_without_warnings
from bitsea.instruments import superfloat as bio_float
from datetime import datetime, timedelta
from bitsea.basins.region import Rectangle
from bitsea.instruments.matchup_manager import Matchup_Manager


inputfile        = args.inputfile
var=args.var
OUTDIR = args.outdir
BASEDIR = args.basedir
TL = TimeList.fromfilenames(None, BASEDIR / "PROFILES/","ave*.nc")
deltaT= timedelta(hours=12)
TI = TimeInterval.fromdatetimes(TL.Timelist[0] - deltaT, TL.Timelist[-1] + deltaT)
ALL_PROFILES = bio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


class ncreader():
    def __init__(self, filename):
        ncIN = NC.Dataset(filename,"r")
        self.nVAR = ncIN.dimensions['var']
        self.nSUB = ncIN.dimensions['sub']
        self.nDEPTH  = ncIN.dimensions['depth']
        self.SUBLIST = ncIN.sublist.split(",")
        self.LAYERLIST=ncIN.layerlist.split(",")
        self.VARLIST = ncIN.varlist.split(",")
        self.npoints = np.array(ncIN.variables['npoints'])
        self.bias    = np.array(ncIN.variables['bias'   ])
        self.rmse    = np.array(ncIN.variables['rmse'   ])
        self.ref    = np.array(ncIN.variables['ref'])
        self.mod    = np.array(ncIN.variables['model'])
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
    
    ref1 = DATAfile.plotdata(DATAfile.ref   , var, sub, layer)
    mod1 = DATAfile.plotdata(DATAfile.mod   , var, sub, layer)

    ii=numb1==0
    rmse1[ii] = np.nan
    bias1[ii] = np.nan

    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()

    ax2.bar(times,numb1,width=7, color='0.5', alpha=0.3, align='center',edgecolor="black", linewidth=3)
    ax2.set_ylabel(' n. of BGC-Argo floats', fontsize=20)
    ax2.set_ylim([0,numb1.max() +2])


    ax.plot(times,bias1,'m.-', label='bias')
    ax.plot(times,rmse1,'k.-', label='rmse')
    if longvar == 'Chlorophyll' : 
        ax.set_ylim([-0.4, 0.4])
        ax.set_ylabel('bias, rmse mg/m$^3$', fontsize=20)
    if longvar == 'PhytoC' :
        ax.set_ylabel('bias, rmse mg/m$^3$', fontsize=20)
    if longvar == 'POC' :
        ax.set_ylabel('bias, rmse mg/m$^3$', fontsize=20)

    if longvar == 'Nitrate'     : 
        ax.set_ylim([-5, 5])
        ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=20)
    if longvar == 'Oxygen'      :
        ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=20)
        
    ax.legend(loc=2)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(16) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(16)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label1.set_fontsize(16)
    for l in ax2.get_yticklabels() : l.set_fontsize(16)

    ii = np.zeros((len(times),) , bool)
    for k,t in enumerate(times) : ii[k] = timeinterval.contains(t)
    biasm = nanmean_without_warnings(bias1[ii])
    rmsem = nanmean_without_warnings(rmse1[ii])
    refm  = nanmean_without_warnings(ref1[ii])
    modm  = nanmean_without_warnings(mod1[ii])
    return fig, biasm, rmsem, ax, ax2, refm, modm

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLONGNAMES={'P_l':'Chlorophyll', 'N3n':'Nitrate', 'O2o':'Oxygen','P_c':'PhytoC', 'POC':'POC'}

SUBLIST = OGS.NRT3.basin_list
nSub = len(SUBLIST)
nLayers = len(LAYERLIST)

ti_restrict = TI

column_names=[layer.string() for layer in LAYERLIST]
row_names   =[sub.name for sub in SUBLIST]


BIAS = np.zeros((nSub,nLayers),np.float32)
RMSE = np.zeros((nSub,nLayers),np.float32)
REF = np.zeros((nSub,nLayers),np.float32)
MOD = np.zeros((nSub,nLayers),np.float32)
for isub, sub in enumerate(SUBLIST):
    for ilayer, layer in enumerate(LAYERLIST):
        outfile = OUTDIR / f"{var}.{sub.name}.{layer.longname()}.png"
        print (outfile,flush=True)
        fig,bias,rmse,ax,ax2, ref, mod  = single_plot(VARLONGNAMES[var],var,sub.name,layer.string(), ti_restrict)
        BIAS[isub,ilayer] = bias
        RMSE[isub,ilayer] = rmse
        REF[isub,ilayer] = ref
        MOD[isub,ilayer] = mod
        title = "%s %s %s " %(VARLONGNAMES[var], sub.extended_name, layer.string())
        fig.suptitle(title, fontsize=20)
        fig.savefig(outfile)
        pl.close(fig)

writetable(OUTDIR /  (var + '_BIAS.txt'),BIAS,row_names, column_names)
writetable(OUTDIR /  (var + '_RMSE.txt'),RMSE,row_names, column_names)
writetable(OUTDIR /  (var + '_REF.txt' ),REF,row_names , column_names)
writetable(OUTDIR /  (var + '_MOD.txt' ),MOD,row_names , column_names)

