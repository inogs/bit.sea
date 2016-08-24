from nrt3_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import numpy as np
from commons.layer import Layer
from basins import V2 as OGS


TI = TimeInterval("20150608","20160101","%Y%m%d")
ARCHIVEDIR="/pico/scratch/userexternal/gbolzon0/NRT/V4/NRT3_outputs/"
V4_data = timelistcontainer(TI,ARCHIVEDIR)

TI = TimeInterval("20160401","20160805","%Y%m%d")
ARCHIVEDIR="/pico/scratch/userexternal/gbolzon0/NRT/V2C/NRT3_outputs/"
V2C_data = timelistcontainer(TI,ARCHIVEDIR)





#times=DL.getTimeList("20150601-12:00:00", "20160701-12:00:00", "days=7")
#nFrames = len(times)



# def single_plot(var):
#     BIAS=np.random.randn(nFrames)
#     RMSE=np.random.randn(nFrames)
#     nPOINTS = np.random.random_integers(0,8,nFrames)
#     ii=nPOINTS==0
#     RMSE[ii] = np.nan
#     BIAS[ii] = np.nan
#
#     pl.close('all')
#     fig, ax = pl.subplots(figsize=(16,4))
#     ax2 = ax.twinx()
#
#     ax2.bar(times,nPOINTS,width=7, color='g', alpha=0.3, label='n points',align='center') #, label='n points'
#     ax2.set_ylabel(' # Points')
#     ax2.set_ylim([0,nPOINTS.max()+2])
#
#     ax.plot(times,BIAS,'r.-' , label='bias')
#     ax.plot(times,RMSE,'b.-', label='rmse')
#     ax.set_ylabel('bias, rmse mg/m$^3$')
#     ax.legend(loc=2)
#     ax2.legend(loc=1)
#     ax.set_title(var)
#
#     return fig

def single_plot(longvar, var, sub, layer ):
    varV4 = var
    if var == 'P_l': varV4 = 'P_i'

    times0, bias0 = V4_data.plotdata(V4_data.bias,   varV4,sub, layer)
    _     , rmse0 = V4_data.plotdata(V4_data.rmse,   varV4,sub, layer)
    _     , numb0 = V4_data.plotdata(V4_data.number, varV4,sub, layer)

    times1, bias1 = V2C_data.plotdata(V2C_data.bias,   var,sub, layer)
    _     , rmse1 = V2C_data.plotdata(V2C_data.rmse,   var,sub, layer)
    _     , numb1 = V2C_data.plotdata(V2C_data.number, var,sub, layer)
    
    ii=numb0==0
    rmse0[ii] = np.nan
    bias0[ii] = np.nan
    ii=numb1==0
    rmse1[ii] = np.nan
    bias1[ii] = np.nan

    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()

    ax2.bar(times0,numb0,width=7, color='g', alpha=0.3, label='n points',align='center')
    ax2.bar(times1,numb1,width=7, color='g', alpha=0.3, align='center')
    ax2.set_ylabel(' # Points')
    ax2.set_ylim([0,max(numb0.max(),numb1.max()) +2])

    ax.plot(times0,bias0,'r.-', label='bias V1')
    ax.plot(times1,bias1,'m.-', label='bias V2')
    ax.plot(times0,rmse0,'b.-', label='rmse V1')
    ax.plot(times1,rmse1,'k.-', label='rmse V2')

    ax.set_ylabel('bias, rmse mg/m$^3$')
    ax.legend(loc=2)
    ax2.legend(loc=1)
    ax.set_title(longvar)

    return fig

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
VARLONGNAMES=['Chlorophyll','Nitrate','Oxygen']
SUBLIST = OGS.NRT3
outputdir='IMG/'

for ivar, var in enumerate(VARLIST):
    for isub, sub in enumerate(SUBLIST):
        for layer in LAYERLIST:
            outfile = "%s%s.%s.%s.png" % (outputdir,var,sub.name,layer.longname())
            print outfile
            fig = single_plot(VARLONGNAMES[ivar],var,sub.name,layer.string())
            fig.savefig(outfile)
            pl.close(fig)
            