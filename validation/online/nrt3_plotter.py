#P_l.tyr.0-100m.png
#N3n.tyr.0-100m.png
#O2o.tyr.0-100m.png


#per ogni png un solo asse in funzione del tempo
#con rmse e bias da una parte, npoints dall'altra
import pylab as pl
import commons.genUserDateList as DL
import numpy as np
from commons.layer import Layer
from basins import V2 as OGS

times=DL.getTimeList("20150601-12:00:00", "20160701-12:00:00", "days=7")
nFrames = len(times)



def single_plot(var):
    BIAS=np.random.randn(nFrames)
    RMSE=np.random.randn(nFrames)
    nPOINTS = np.random.random_integers(0,8,nFrames)
    ii=nPOINTS==0
    RMSE[ii] = np.nan
    BIAS[ii] = np.nan
    
    pl.close('all')
    fig, ax = pl.subplots(figsize=(16,4))
    ax2 = ax.twinx()
    
    ax2.bar(times,nPOINTS,width=7, color='g', alpha=0.3, label='n points',align='center') #, label='n points'
    ax2.set_ylabel(' # Points')
    ax2.set_ylim([0,nPOINTS.max()+2])
    
    ax.plot(times,BIAS,'r.-' , label='bias')
    ax.plot(times,RMSE,'b.-', label='rmse')
    ax.set_ylabel('bias, rmse mg/m$^3$')
    ax.legend(loc=2)
    ax2.legend(loc=1)
    ax.set_title(var)

    return fig

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
VARLONGNAMES=['Chlorophyll','Nitrate','Oxygen']
outputdir='IMG/'
for ivar, var in enumerate(VARLIST):
    for isub, sub in enumerate(OGS.NRT3):
        for layer in LAYERLIST:
            outfile = "%s%s.%s.%s.png" % (outputdir,var,sub.name,layer.longname())
            print outfile
            fig = single_plot(VARLONGNAMES[ivar])
            fig.savefig(outfile)
            pl.close(fig)
            