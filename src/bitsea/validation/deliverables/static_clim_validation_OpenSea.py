import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates in output directory two files (model and ref) 
    containing [nSub, nLayers] arrays of climatologies, 
    both for the standard 8 layers and for 14 layers.
    These arrays will be used in the next step to generate the following metrics:
    N1p-LAYER-Y-CLASS4-CLIM
    N1p-LAYER-Y-CLASS4-CLIM_STD
    N1p-PROF-Y-CLASS4-CLIM-CORR-BASIN
    and analogous ones for N3n, N4n, N5s, O2o, ALK, DIC, pH, pCO2
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = '')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--climdir', '-c',
                                type = str,
                                default = "/g100_work/OGS_test2528/Observations/TIME_RAW_DATA/STATIC/MedBGCins",
                                required = True,
                                help = "Directory with Clim_Annual_Ref inside containing NetCDF files")
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--starttime','-s',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-e',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')   
    return parser.parse_args()

args = argument()

import numpy as np
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.timeseries.plot import Hovmoeller_matrix
from bitsea.timeseries.plot import read_pickle_file
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer
from bitsea.basins import V2 as basV2
from bitsea.commons.utils import addsep
from bitsea.matchup.statistics import matchup
from bitsea.commons.utils import writetable
from bitsea.commons import timerequestors
import xarray as xr

#List consistent with /g100_work/OGS_test2528/Observations/TIME_RAW_DATA/STATIC/MedBGCins/Clim_Annual_Ref
LayerList = [Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]

INPUTDIR=addsep(args.inputdir)
OUTDIR = addsep(args.outdir)
CLIMDIR = addsep(args.climdir) + "Clim_Annual_Ref/"

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")

TheMask = Mask.from_file(args.maskfile)
jpk,jpj,jpi = TheMask.shape


TL = TimeList.fromfilenames(TI, INPUTDIR, "ave*nc")
Req=timerequestors.Generic_req(TI)
ind,ww=TL.select(Req) 

def Layers_Mean(Pres,Values,LayerList):
    '''
    Performs mean of profile along layers.

    Arguments:
    * Pres      * numpy array of pressures
    * Values    * numpy array of concetrations
    * LayerList * list of layer objects
    Returns :
    * MEAN_LAY * [nLayers] numpy array, mean of each layer
    '''
    MEAN_LAY = np.zeros(len(LayerList), np.float32)

    for ilayer, layer in enumerate(LayerList):
        ii = (Pres>=layer.top) & (Pres<layer.bottom)
        if (ii.sum()> 1 ) :
            local_profile = np.array(Values[ii])
            if not np.isnan(local_profile).any():
                MEAN_LAY[ilayer] = np.mean(local_profile)
    return MEAN_LAY


VARLIST=['N1p','N3n','O2o','ALK','DIC','pH','N4n','N5s','pCO2']
SUBlist = basV2.Pred.basin_list[:]
nLayers = len(LayerList)

METRICvar = {'N1p':'PHO',
             'N3n':'NIT',
             'O2o':'DO',
             'ALK':'ALK',
             'DIC':'DIC',
             'pH':'pH_t',
             'N4n':'NH4',
             'N5s':'SiO2',
             'pCO2':'pCO2'}

# remove Atlantic buffer from the list
SUBlist.remove(SUBlist[-1])

rows_names  =[layer.string() for layer in LayerList]
column_names=['bias','rmse','corr']
column_names_STD=['bias','rmse','corr','mod_MEAN','ref_MEAN','mod_STD','ref_STD']
for ivar, var in enumerate(VARLIST):
    filename = INPUTDIR + var + ".pkl"
    TIMESERIES_complete,TL_complete=read_pickle_file(filename)
    ind,ww=TL_complete.select(Req) 
    TIMESERIES=TIMESERIES_complete[ind,:]
    climfile= CLIMDIR + var + "_clim_metrics.nc"
    with xr.open_dataset(climfile) as ds:
        CLIM_REF_static = ds['mean'].values

    nSub = len(SUBlist)
    CLIM_MODEL = np.zeros((nSub, nLayers))
    for iSub, sub in enumerate(SUBlist):
        Mean_profiles,_,_ = Hovmoeller_matrix(TIMESERIES,TL, np.arange(jpk), iSub, icoast=1, istat=0)
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        CLIM_MODEL[iSub,:] = Layers_Mean(TheMask.zlevels, mean_profile,LayerList)
    np.save(OUTDIR + var + "_ref_clim", CLIM_REF_static)
    np.save(OUTDIR + var + "_mod_clim", CLIM_MODEL)
    STATS = np.zeros((nLayers,3),np.float32)*np.nan
    STATS_STD = np.zeros((nLayers,7),np.float32)*np.nan
    for ilayer, layer in enumerate(LayerList):
        refsubs = CLIM_REF_static[:,ilayer]
        modsubs =      CLIM_MODEL[:,ilayer]
        bad = np.isnan(refsubs) | np.isnan(modsubs)
        good = ~bad

        m = matchup(modsubs[good], refsubs[good])

        STATS[ilayer,0] = m.bias()
        STATS[ilayer,1] = m.RMSE()
        STATS[ilayer,2] = m.correlation()
 
        STATS_STD[ilayer,0] = m.bias()
        STATS_STD[ilayer,1] = m.RMSE()
        STATS_STD[ilayer,2] = m.correlation()
        STATS_STD[ilayer,3] = np.mean(modsubs[good])
        STATS_STD[ilayer,4] = np.mean(refsubs[good])
        STATS_STD[ilayer,5] = np.std(modsubs[good])
        STATS_STD[ilayer,6] = np.std(refsubs[good])
    writetable(OUTDIR + var + "-LAYER-Y-CLASS4-CLIM.txt", STATS,rows_names,column_names)
    writetable(OUTDIR + var + "-LAYER-Y-CLASS4-CLIM_STD.txt", STATS_STD,rows_names,column_names_STD)

# PROF statistics
def clim_prof_Qc(CLIMDIR, var):
    """
    Applies quality check flag to climatology
    profiles used in plot (14 layers)
    """
    climfile= CLIMDIR + var + "_clim_plot.nc"
    with xr.open_dataset(climfile) as ds:
        clim_tmp = ds['mean'].values
    climfile= CLIMDIR + var + "_clim_metrics.nc"
    with xr.open_dataset(climfile) as ds:
        included = ds['included_profs_flag'].values
    qc = included.astype(np.float64)
    qc[qc == 0] = np.nan
    CLIM_REF_static = clim_tmp*qc[:,None]
    return CLIM_REF_static

PresDOWN=np.array([25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList_2=[]
top = 0
for bottom in PresDOWN:
    LayerList_2.append(Layer(top, bottom))
    top = bottom
nLayers = len(LayerList_2)
LayerList_3=LayerList_2[:7]
nLayers3 = len(LayerList_3)

rows_names=[sub.name for sub in SUBlist]
column_names = ['correlation']


for var in VARLIST:
    if ( var == "pCO2" ): # Integral only up to 200m. It works because pCO2 is the last var of the VARLIST 
       LayerList_2 = LayerList_3
       nLayers = nLayers3

    filename = INPUTDIR + var + ".pkl"
    TIMESERIES_complete,TL_complete=read_pickle_file(filename)
    ind,ww=TL_complete.select(Req)
    TIMESERIES=TIMESERIES_complete[ind,:]
    CLIM_REF_static = clim_prof_Qc(CLIMDIR,var)

    if ( var == "pCO2" ):
        CLIM_REF_static = CLIM_REF_static[:,:7]


    CLIM_MODEL = np.zeros((nSub, nLayers))
    for iSub, sub in enumerate(SUBlist):
        Mean_profiles,_,_ = Hovmoeller_matrix(TIMESERIES,TL, np.arange(jpk), iSub, icoast=1, istat=0)
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        CLIM_MODEL[iSub,:] = Layers_Mean(TheMask.zlevels, mean_profile,LayerList_2)
    np.save(OUTDIR + var + "ref_clim14", CLIM_REF_static)
    np.save(OUTDIR + var + "mod_clim14", CLIM_MODEL)
    STATS = np.zeros((nSub,3),np.float32)*np.nan
    for iSub, sub in enumerate(SUBlist):
        refsubs = CLIM_REF_static[iSub,:]
        modsubs =      CLIM_MODEL[iSub,:]
        bad = np.isnan(refsubs) | np.isnan(modsubs)
        good = ~bad
        ngoodlayers=good.sum()
        if ngoodlayers>0:
            m = matchup(modsubs[good], refsubs[good])
            STATS[iSub,0] = m.bias()
            STATS[iSub,1] = m.RMSE()
            STATS[iSub,2] = m.correlation()
    writetable(OUTDIR + var + "-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt", STATS, rows_names, ['bias','rmse','corr'])

