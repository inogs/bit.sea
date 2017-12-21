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
    parser.add_argument(   '--validation_dir','-v',
                                type = str,
                                required = True,
                                help = 'wrkdir/2/POSTPROC/AVE_FREQ_1/online_validation')
    parser.add_argument(   '--previous_archive','-p',
                                type = str,
                                required = True,
                                help = 'previous chain archive directory, taken from static-data')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()


from biofloats_ms_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import numpy as np
from commons.layer import Layer
from basins import V2 as OGS
from commons.utils import addsep
from dateutil.relativedelta import relativedelta
from datetime import datetime

ARCHIVEDIR       = addsep(args.archivedir)
ARCHIVE_PREV     = addsep(args.previous_archive)
VALID_WRKDIR     = addsep(args.validation_dir)
OUTFIG_DIR       = addsep(args.outdir)

Graphic_DeltaT = relativedelta(months=18)
datestart = datetime.strptime(args.date,'%Y%m%d') -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")

#TI_V2 = TimeInterval("20160401","20171206","%Y%m%d")
TI_V2 = TimeInterval(timestart,"20171206","%Y%m%d")
TI_V3 = TimeInterval('20171114', args.date,'%Y%m%d')

#V4_data  = timelistcontainer(TI_V1,ARCHIVE_PREV,postfix_dir="")
#V2C_data = timelistcontainer(TI_V2,ARCHIVEDIR  ,postfix_dir="POSTPROC/AVE_FREQ_1/validation/biofloats_ms/")
#V2C_data = timelistcontainer(TI_V2,ARCHIVEDIR,postfix_dir="dir1/dir2/")
V2_data = timelistcontainer(TI_V2,ARCHIVE_PREV,postfix_dir="")
V3_data = timelistcontainer(TI_V3,ARCHIVEDIR,postfix_dir="POSTPROC/AVE_FREQ_1/validation/biofloats/")
V3_data.append_dir(VALID_WRKDIR)


def single_plot(longvar, var, sub, layer ):
    varV2 = var
#    if var == 'P_l': varV4 = 'P_i'

    times0, bias0 = V2_data.plotdata(V2_data.bias,   varV2,sub, layer)
    _     , rmse0 = V2_data.plotdata(V2_data.rmse,   varV2,sub, layer)
    _     , numb0 = V2_data.plotdata(V2_data.number, varV2,sub, layer)

    times1, bias1 = V3_data.plotdata(V3_data.bias,   var,sub, layer)
    _     , rmse1 = V3_data.plotdata(V3_data.rmse,   var,sub, layer)
    _     , numb1 = V3_data.plotdata(V3_data.number, var,sub, layer)
    
    ii=numb0==0
    rmse0[ii] = np.nan
    bias0[ii] = np.nan
    ii=numb1==0
    rmse1[ii] = np.nan
    bias1[ii] = np.nan

    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()

    ax2.bar(times0,numb0,width=7, color='0.3', alpha=0.3, align='center', edgecolor="k")
    ax2.bar(times1,numb1,width=7, color='0.3', alpha=0.3, align='center', edgecolor="k")
    ax2.set_ylabel(' n. of BGC-Argo floats', fontsize=20)
    ax2.set_ylim([0,max(numb0.max(),numb1.max()) +2])
#    ax2.set_ylim([0,numb1.max() +2])

    ax.plot(times0,bias0,'r.-', label='bias V2')
    ax.plot(times1,bias1,'m.-', label='bias V3')
    ax.plot(times0,rmse0,'b.-', label='rmse V2')
    ax.plot(times1,rmse1,'k.-', label='rmse V3')

#    ax.set_ylabel('bias, rmse mg/m$^3$')

    if longvar == 'Chlorophyll' :
        ax.set_ylim([-0.4, 0.4])
        ax.set_ylabel('bias, rmse mg/m$^3$', fontsize=20)
#       ax.ticklabel_format(axis='both', fontsize=20)

    if longvar == 'Nitrate'     :
        ax.set_ylim([-5, 5])
        ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=20)
    if longvar == 'Oxygen'      :
        ax.set_ylabel('bias, rmse mmol/m$^3$', fontsize=20)

    ax.legend(loc=2)
    ax2.legend(loc=1)
    ax.set_title(longvar)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for l in ax2.get_yticklabels() : l.set_fontsize(16)


#    ii = np.zeros((len(times),) , np.bool)
#    for k,t in enumerate(times) : ii[k] = timeinterval.contains(t)
#    biasm = np.nanmean(bias1[ii])
#    rmsem = np.nanmean(rmse1[ii])
#    return fig, biasm, rmsem, ax, ax2
#    return fig, ax, ax2
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

if __name__ == '__main__':
    from commons.time_interval import TimeInterval
    TI_V2 = TimeInterval("20160412","20170502","%Y%m%d")
    A_DIR="/pico/scratch/userexternal/lfeudale/ANALYSIS/NRT3_outputs/"
    V_WRKDIR = "/pico/scratch/userexternal/lfeudale/ANALYSIS/NRT3_outputs/" 
    FIG_DIR = "/pico/scratch/userexternal/lfeudale/ANALYSIS/NRT3_figs/"
#    run biofloats_ms_plotter.py -d '20170502' -a A_DIR -v V_WRKDIR -o FIG_DIR	

