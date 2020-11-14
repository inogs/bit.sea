import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates in output directory two files ( model and ref) 
    containing [nSub, nLayers] arrays of climatologies.
    These arrays will be used in the next step to generate the following metrics:

    PHO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    NIT-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
     DO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    ALK-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    DIC-LAYER-Y-CLASS4-CLIM-BIAS/RMSD

    ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN
    DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN
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
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from timeseries.plot import Hovmoeller_matrix
from timeseries.plot import read_pickle_file
from commons.mask import Mask
from commons.layer import Layer
from basins import V2 as basV2
from static import climatology
from commons.utils import addsep
from matchup.statistics import matchup
from commons.utils import writetable

LayerList = [Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]

INPUTDIR=addsep(args.inputdir)
OUTDIR = addsep(args.outdir)
TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")

TheMask= Mask(args.maskfile, loadtmask=False)
jpk,jpj,jpi = TheMask.shape
z = -TheMask.zlevels

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

TL = TimeList.fromfilenames(TI, INPUTDIR, "ave*nc")

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
        ii = (Pres>=layer.top) & (Pres<=layer.bottom)
        if (ii.sum()> 1 ) :
            local_profile = Values[ii]
            MEAN_LAY[ilayer] = np.mean(local_profile)
    return MEAN_LAY

# BFMv2:
#VARLIST=['N1p','N3n','O2o','Ac','DIC','pH']
# BFMv5:
VARLIST=['N1p','N3n','O2o','ALK','DIC','pH','pCO2','N4n','N5s']
SUBlist = basV2.Pred.basin_list
#SUBlist2 = basV2.Pred2.basin_list
#nSub = len(SUBlist)
nLayers = len(LayerList)
METRICvar = {'N1p':'PHO',
             'N3n':'NIT',
             'O2o':'DO',
             'Ac':'ALK',
             'ALK':'ALK',
             'DIC':'DIC',
             'pH':'pH_t',
            'pCO2':'pCO2',
             'N4n':'NH4',
             'N5s':'SiO2'}


# Remove Altantic Buffer from the list:
SUBlist.remove(SUBlist[-1])

rows_names  =[layer.string() for layer in LayerList]
column_names=['bias','rmse','corr']
column_names_STD=['bias','rmse','corr','mod_MEAN','ref_MEAN','mod_STD','ref_STD']
for ivar, var in enumerate(VARLIST):
#  if (ivar == 4) : #3 7 
    filename = INPUTDIR + var + ".pkl"
    TIMESERIES,TL=read_pickle_file(filename)
    print METRICvar[var] + "-LAYER-Y-CLASS4-CLIM-BIAS,RMSD"
#    if ( var in ["N1p","N3n","N5s","O2o","O3c","O3h"] ):
#    CLIM_REF_static,_ = climatology.get_climatology(var,SUBlist, LayerList, basin_expand=True)
#    else:
#        CLIM_REF_static,_ = climatology.get_climatology(var,SUBlist, LayerList)

    
#    if ( var == 'N1p' ): 
#        SUBlist = basV2.Pred2.basin_list
#        print len(SUBlist)
#    else: SUBlist = basV2.Pred.basin_list
    CLIM_REF_static,_ = climatology.get_climatology(var,SUBlist, LayerList, basin_expand=True,QC=True)
    nSub = len(SUBlist)
    CLIM_MODEL = np.zeros((nSub, nLayers))
    for iSub, sub in enumerate(SUBlist):
        print sub.name
        Mean_profiles,_,_ = Hovmoeller_matrix(TIMESERIES,TL, np.arange(jpk), iSub, icoast=1, istat=0)
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        CLIM_MODEL[iSub,:] = Layers_Mean(TheMask.zlevels, mean_profile,LayerList)
    np.save(OUTDIR + var + "_ref_clim", CLIM_REF_static)
    np.save(OUTDIR + var + "_mod_clim", CLIM_MODEL)
    print "CLIM_REF_static"
    print CLIM_REF_static
    print "CLIM_MODEL"
    print CLIM_MODEL
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

# N1p e N3n in table 4.6
# O2o in table 4.10
# Alk, dic, pH, pCO2 in table 4.13



PresDOWN=np.array([25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList_2=[]
top = 0
for bottom in PresDOWN:
    LayerList_2.append(Layer(top, bottom))
    top = bottom
nLayers = len(LayerList_2)
LayerList_3=LayerList_2[:7]
print "LayerList_3 = ", LayerList_3
nLayers3 = len(LayerList_3)


rows_names=[sub.name for sub in SUBlist]
column_names = ['correlation']

for var in VARLIST:
    if ( var == "pCO2" ): # Integral only up to 200m. It works because pCO2 is the last var of the VARLIST 
       LayerList_2 = LayerList_3
       nLayers = nLayers3

    filename = INPUTDIR + var + ".pkl"
    TIMESERIES,TL=read_pickle_file(filename)
    print METRICvar[var] + "-PROF-Y-CLASS4-CLIM-CORR-BASIN"
#    if ( var in ["N1p","N3n","N5s","O2o","O3c","O3h"] ):
    CLIM_REF_static,_ = climatology.get_climatology(var,SUBlist, LayerList_2, basin_expand=True,QC=True)
#    else:
#        CLIM_REF_static,_ = climatology.get_climatology(var,SUBlist, LayerList_2)

    CLIM_MODEL = np.zeros((nSub, nLayers))
    for iSub, sub in enumerate(SUBlist):
        Mean_profiles,_,_ = Hovmoeller_matrix(TIMESERIES,TL, np.arange(jpk), iSub, icoast=1, istat=0)
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        CLIM_MODEL[iSub,:] = Layers_Mean(TheMask.zlevels, mean_profile,LayerList_2)
    np.save(OUTDIR + var + "ref_clim14", CLIM_REF_static)
    np.save(OUTDIR + var + "mod_clim14", CLIM_MODEL)
#    CORR = np.zeros((nSub,1),np.float32)*np.nan
    STATS = np.zeros((nSub,3),np.float32)*np.nan
    for iSub, sub in enumerate(SUBlist):
        refsubs = CLIM_REF_static[iSub,:]
        modsubs =      CLIM_MODEL[iSub,:]
        bad = np.isnan(refsubs) | np.isnan(modsubs)
        good = ~bad
        ngoodlayers=good.sum()
        if ngoodlayers>0:
            m = matchup(modsubs[good], refsubs[good])
#            CORR[iSub,0] = m.correlation()
	    STATS[iSub,0] = m.bias()
	    STATS[iSub,1] = m.RMSE()
	    STATS[iSub,2] = m.correlation()
#    writetable(OUTDIR + var + "-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt", CORR, rows_names,column_names)
    writetable(OUTDIR + var + "-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt", STATS, rows_names, ['bias','rmse','corr'])

# Table 4.7 Correlazione N1p, N3n per certi subbasins
# Table 4.11 Correlazione O2o per certi subbasins
# Table 4.14 Bias,Rms corr per V2.Pred subs ALK, DIC, pH, pCO2
