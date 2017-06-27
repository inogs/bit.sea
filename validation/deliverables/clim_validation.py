import numpy as np
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from timeseries.plot import Hovmoeller_matrix
from commons.mask import Mask
from commons.layer import Layer
from basins import V2 as basV2
from static import climatology

PresDOWN=np.array([10,30,60,100,150,300,600,1000])
LayerList=[]
top = 0
for bottom in PresDOWN:
    LayerList.append(Layer(top, bottom))
    top = bottom
    

MODDIR="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
TI = TimeInterval("20140101","20180101","%Y%m%d")
maskfile ="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc"
OUTDIR = "static_clim/"

TheMask= Mask(maskfile, loadtmask=False)
jpk,jpj,jpi = TheMask.shape
z = -TheMask.zlevels

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")

def Layers_Mean(Pres,Values):
    MEAN_LAY = np.zeros(len(LayerList), np.float32)

    for ilayer, layer in enumerate(LayerList):
        ii = (Pres>=layer.top) & (Pres<=layer.bottom)
        if (ii.sum()> 1 ) :
            local_pres = Pres[ii]
            local_profile = Values[ii]
            MEAN_LAY[ilayer] = np.mean(local_profile)
    return MEAN_LAY

VARLIST=['N1p','N3n','O2o','Ac','DIC']
nSub = len(basV2.P.basin_list)
nLayers = len(LayerList)


for ivar, var in enumerate(VARLIST):
    print var
    filename = OUTDIR + var + "ref_clim"
    CLIM_REF_static = climatology.get_climatology(var,basV2.P.basin_list, LayerList)
    
    CLIM_MODEL = np.zeros((nSub, nLayers))
    for iSub, sub in enumerate(basV2.P):
        Mean_profiles,_,_ = Hovmoeller_matrix(TL.Timelist,TL.filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        CLIM_MODEL[iSub,:] = Layers_Mean(TheMask.zlevels, mean_profile)
    np.save(OUTDIR + var + "ref_clim", CLIM_REF_static)
    np.save(OUTDIR + var + "mod_clim", CLIM_MODEL)

