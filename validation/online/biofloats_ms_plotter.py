import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates time series png files for mdeeaf web site.
    Each file has 3 subplots: bias, rms and number of points.
    At each run time becomes longer.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--date','-d',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')

    parser.add_argument(   '--archivedir','-a',
                                type = str,
                                required = True,
                                help = 'chain archive directory')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')
from biofloats_ms_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
import matplotlib.pyplot as pl
import numpy as np
from commons.layer import Layer
from basins import V2 as OGS
from commons.utils import addsep
from dateutil.relativedelta import relativedelta
from datetime import datetime
from matplotlib.ticker import MaxNLocator

ARCHIVEDIR       = addsep(args.archivedir)
OUTFIG_DIR       = addsep(args.outdir)

Graphic_DeltaT = relativedelta(months=18)
datestart = datetime.strptime(args.date,'%Y%m%d') -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")


TI = TimeInterval('20171114', args.date,'%Y%m%d')

prefix="BioFloat_Weekly_validation_"
data = timelistcontainer(TI,ARCHIVEDIR,"BioFloat_Weekly_validation_*", prefix=prefix)
#data = timelistcontainer(TI,ARCHIVEDIR,postfix_dir="POSTPROC/AVE_FREQ_1/validation/biofloats/")


def single_plot(longvar, var, sub, layer ):
    varV2 = var

    times1, bias1 = data.plotdata(data.bias,   var,sub, layer)
    _     , rmse1 = data.plotdata(data.rmse,   var,sub, layer)
    _     , numb1 = data.plotdata(data.number, var,sub, layer)
    

    ii=numb1==0
    rmse1[ii] = np.nan
    bias1[ii] = np.nan

    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()


    ax2.bar(times1,numb1,width=7, color='0.9', align='center', edgecolor="k")
    ax2.set_ylabel(' n. of BGC-Argo floats', fontsize=16)
    ax2.set_ylim([0,numb1.max() +2])
#    ax2.set_ylim([0,numb1.max() +2])


    ax.plot(times1,bias1,'m.-', label='bias')
    ax.plot(times1,rmse1,'k.-', label='rmse')

#    ax.set_ylabel('bias, rmse mg/m$^3$')

    if longvar == 'Chlorophyll' :
        ax.set_ylim([-0.4, 0.4])
        ax.set_ylabel('bias, rmse mg/m$^3$', fontsize=16)
#       ax.ticklabel_format(axis='both', fontsize=20)

    if longvar == 'Nitrate'     :
        ax.set_ylim([-5, 5])
        ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=16)
    if longvar == 'Oxygen'      :
        ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=16)

    ax.legend(loc=2)
    ax.set_title(longvar, fontsize=20)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    for l in ax2.get_yticklabels() : l.set_fontsize(16)


    ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
    ax.patch.set_visible(False) # hide the 'canvas'
    return fig

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
VARLONGNAMES=['Chlorophyll','Nitrate','Oxygen']
SUBLIST = OGS.NRT3


for ivar, var in enumerate(VARLIST):
    for isub, sub in enumerate(SUBLIST):
        for layer in LAYERLIST:
            outfile = "%s%s.%s.%s.png" % (OUTFIG_DIR,var,sub.name,layer.longname())
            print outfile
            fig = single_plot(VARLONGNAMES[ivar],var,sub.name,layer.string())
            fig.savefig(outfile)
            pl.close(fig)

