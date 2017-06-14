from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from commons.submask import SubMask
from commons.layer import Layer
from timeseries.plot import Hovmoeller_matrix
import numpy as np
from basins import V2
import pylab as pl
from commons.utils import getcolor, addsep

def weighted_mean(values,weights):
    if weights.sum() == 0 :
        return np.nan
    else:
        return (m*weights).sum()/(weights.sum())

EAST_SUBBASIN_LIST=[V2.lev1,V2.lev2,V2.lev3,V2.lev4,V2.aeg]
WEST_SUBBASIN_LIST=[V2.alb, V2.nwm,V2.tyr1,V2.tyr2,V2.swm1, V2.swm2]
MID__SUBBASIN_LIST=[V2.adr1, V2.adr2, V2.ion1,V2.ion2,V2.ion3]

region_list = [WEST_SUBBASIN_LIST,MID__SUBBASIN_LIST,EAST_SUBBASIN_LIST]
region_names= ['West','Central','East']


LAYERLIST = [Layer(0,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300),Layer(300,600)]
nLayers = len(LAYERLIST)

def figure_generator(subbasin_list,var,region_name):
    fig,axes=pl.subplots(nLayers,1)
    fig.set_size_inches(16,10)
    fig.set_dpi(150)
    axes=axes.ravel()
    x_shift = -0.025
    for ilayer,layer in enumerate(LAYERLIST):
        ax=axes[ilayer]
        ax.set_ylabel(layer.string())
        BBox = ax.get_position()
        rect = [BBox.xmin+x_shift, BBox.ymin, BBox.width, BBox.height]
        ax.set_position(rect)
        ax.locator_params(axis='y', nbins=4)
    for ax in axes[:nLayers-1]:
        ax.set_xticklabels([])
    title = "%s  %s Mediterranean Sea"  % (var,region_name)
    fig.suptitle(title)
    return fig, axes



TheMask=Mask('/Users/gbolzon/Documents/workspace/ogs_bounday_conditions/masks/meshmask_872.nc')
jpk,jpj,jpi=TheMask.shape
area = TheMask.area
e3t = TheMask.dz
Volume=np.zeros((jpk,jpj,jpi),dtype=np.double)
for i in range(jpk):
    Volume[i,:,:] = area*e3t[i]


nSub    = len(V2.P.basin_list)
LAYER_VOLUME = np.zeros((jpk,nSub), np.float32)
for iSub, sub in enumerate(V2.P.basin_list):
    submask=SubMask(sub,maskobject = TheMask).mask
    for k in range(jpk):
        V=Volume[k,:,:]
        m=submask[k,:,:]
        LAYER_VOLUME[k,iSub] = V[m].sum()


INPUTDIR="/Users/gbolzon/Downloads/STAT_PROFILES"
OUTDIR = addsep("/Users/gbolzon/Documents/workspace/bit.sea/validation/deliverables/out")
TI = TimeInterval("201401","201601","%Y%m")
TL =TimeList.fromfilenames(TI, INPUTDIR, "ave*nc")


VARLIST=['N1p','N3n','O2o','N5s','ppn','P_l','pCO2', 'Ac', 'B1c', 'P_c','R2c','pH','DIC']






for var in VARLIST[:3]:
    for iregion, SUBBASIN_LIST in enumerate(region_list):
        region = region_names[iregion]
        outfile = OUTDIR + var + "." + region + ".png"
        print outfile
        fig, axes = figure_generator(SUBBASIN_LIST, var, region)
        for iSub, sub in enumerate(SUBBASIN_LIST):
            overall_isub = V2.P.basin_list.index(sub)
            Mean_profiles,_,_ = Hovmoeller_matrix(TL.Timelist,TL.filelist, var, overall_isub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
            for ilayer,layer in enumerate(LAYERLIST):
                ax=axes[ilayer]
                ii = (TheMask.zlevels>=layer.top) & (TheMask.zlevels<=layer.bottom)
                weights = LAYER_VOLUME[ii,overall_isub]
                y = np.zeros((TL.nTimes,))
                for it in range(TL.nTimes):
                    m = Mean_profiles[ii,it]
                    y[it] = weighted_mean(m, weights)
                color=getcolor(len(SUBBASIN_LIST), iSub)
                ax.plot(TL.Timelist,y,color=color,label=sub.name)
                
        ax=axes[0]        
        ax.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0,labelspacing=0.25, handletextpad=0,borderpad=1.6)
        leg = ax.get_legend()
        ltext  = leg.get_texts()
        pl.setp(ltext,fontsize=11)
        fig.savefig(outfile)
        pl.close(fig)


    
    
