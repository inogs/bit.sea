import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png files for fig4.11 and fig4.17
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

import matplotlib
matplotlib.use('Agg')

from commons.layer import Layer
import basins.V2 as basV2
from static.climatology import get_climatology_coast
import figure_generator
import figure_generator_extended as fg2
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from timeseries.plot import Hovmoeller_matrix
from timeseries.plot import read_pickle_file, read_basic_info
from commons import timerequestors
import numpy as np
from commons.mask import Mask
from commons.submask import SubMask
import matplotlib.pyplot as pl
from commons.utils import addsep
from commons import genUserDateList as DL

IDrun='new_dset'
OUTDIR=addsep(args.outdir)
MODDIR=addsep(args.inputdir)

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
Req=timerequestors.Generic_req(TI)

TheMask = Mask.from_file(args.maskfile)
jpk,jpj,jpi = TheMask.shape
z = -TheMask.zlevels

#List consistent with /g100_work/OGS_test2528/Observations/TIME_RAW_DATA/STATIC/MedBGCins/Clim_Annual_Ref
PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

# In case of quid production, use the row below
#TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")
weekly=DL.getTimeList(TI.start_time, TI.end_time, days=7)
TL=TimeList(weekly)

SUBLIST = basV2.P.basin_list

N3n_clim, N3n_std = get_climatology_coast('N3n', SUBLIST, LayerList,TheMask,basin_expand=False)
N1p_clim, N1p_std = get_climatology_coast('N1p', SUBLIST, LayerList,TheMask,basin_expand=False)
O2o_clim, O2o_std = get_climatology_coast('O2o', SUBLIST, LayerList,TheMask,basin_expand=False)
N4n_clim, N4n_std = get_climatology_coast('N4n', SUBLIST, LayerList,TheMask,basin_expand=False)
N5s_clim, N5s_std = get_climatology_coast('N5s', SUBLIST, LayerList,TheMask,basin_expand=False)



VARLIST=['P_l','N1p','N3n','O2o','N4n','N5s']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

# reading of input files is limited here ----------
_, COASTLIST, STAT_LIST = read_basic_info('/g100_scratch/userexternal/vdibiagi/EMODnet_2022/forValidation/ave.20190110-12:00:00.stat_profiles.nc')
icoast = COASTLIST.index('coast')
istat =  STAT_LIST.index('Mean')

timeseries_DICT={}
for var in VARLIST:
    filename=MODDIR + var + ".pkl"
    TIMESERIES_complete,TL_complete=read_pickle_file(filename)
    ind,ww=TL_complete.select(Req)
    TIMESERIES=TIMESERIES_complete[ind,:]
    timeseries_DICT[var]=TIMESERIES
#-------------------------------------------------

for iSub, sub in enumerate(basV2.P):
#  if (1==0):
    submask = SubMask(sub, TheMask)
    F = fg2.figure_generator(submask)
    fig, axes = F.gen_structure_1(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_4.11_annual." + sub.name + ".png"
    
    for iTime , tTime in enumerate(TL.Timelist):
        datetimelist= [TL.Timelist[iTime] ]
        MEAN = np.zeros((jpk,), dtype=var_dype )
        for ivar, var in enumerate(VARLIST):
            TIMESERIESvar=timeseries_DICT[var]
            mean_profile = TIMESERIESvar[iTime,iSub,icoast,:,istat]
            #Mean_profiles,_,_ = Hovmoeller_matrix(datetimelist,[filename], var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
            #mean_profile = Mean_profiles[:,0]#.mean(axis=1)
            mean_profile[mean_profile==0]=np.nan
            MEAN[var] = mean_profile
        
        label = 'none'
        color = '0.4'
        fg2.profile_plotter(z,MEAN['P_l'],color, axes[0], None,   label)
        fg2.profile_plotter(z,MEAN['N3n'],color, axes[1], axes[7],label)
        fg2.profile_plotter(z,MEAN['N1p'],color, axes[2], axes[8],label)
        fg2.profile_plotter(z,MEAN['O2o'],color, axes[3], axes[9],label)            
        fg2.profile_plotter(z,MEAN['N5s'],color, axes[4], axes[10],label) 
        fg2.profile_plotter(z,MEAN['N4n'],color, axes[5], axes[11],label)
        MEAN = np.zeros((jpk,), dtype=var_dype )

    for ivar, var in enumerate(VARLIST):
        TIMESERIESvar=timeseries_DICT[var]
        mean_profile = TIMESERIESvar[:,iSub,icoast,:,istat].mean(axis=0)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile    
    fg2.profile_plotter(z,MEAN['P_l'],'k', axes[0], None,   label)
    fg2.profile_plotter(z,MEAN['N3n'],'k', axes[1], axes[7],label)
    fg2.profile_plotter(z,MEAN['N1p'],'k', axes[2], axes[8],label)
    fg2.profile_plotter(z,MEAN['O2o'],'k', axes[3], axes[9],label) 
    fg2.profile_plotter(z,MEAN['N5s'],'k', axes[4], axes[10],label) 
    fg2.profile_plotter(z,MEAN['N4n'],'k', axes[5], axes[11],label)
    
    N3n_clim_mean=N3n_clim[iSub,:]
    N3n_clim_std =N3n_std[ iSub,:]
    N1p_clim_mean=N1p_clim[iSub,:]
    N1p_clim_std =N1p_std[ iSub,:]
    O2o_clim_mean=O2o_clim[iSub,:]
    O2o_clim_std =O2o_std[ iSub,:]
    N5s_clim_mean=N5s_clim[iSub,:]
    N5s_clim_std =N5s_std[ iSub,:]
    N4n_clim_mean=N4n_clim[iSub,:]
    N4n_clim_std =N4n_std[ iSub,:]
    
    fg2.clim_profile_plotter(z_clim,N3n_clim_mean,N3n_clim_std, axes[1], axes[7])
    fg2.clim_profile_plotter(z_clim,N1p_clim_mean,N1p_clim_std, axes[2], axes[8])
    fg2.clim_profile_plotter(z_clim,O2o_clim_mean,O2o_clim_std, axes[3], axes[9])
    fg2.clim_profile_plotter(z_clim,N5s_clim_mean,N5s_clim_std, axes[4], axes[10])
    fg2.clim_profile_plotter(z_clim,N4n_clim_mean,N4n_clim_std, axes[5], axes[11])

    fig.savefig(outfile)
    print (outfile,flush=True)
    pl.close(fig)


# Figures 4.17 (Previously Fig4.19)
#try:
#    Ac__clim, Ac__std = get_climatology_coast('Ac' , SUBLIST, TheMask, LayerList)
#except:
#    print "no adr1 for ALK"
#    Ac__clim, Ac__std = get_climatology_coast('Ac' , SUBLIST2, TheMask, LayerList)
##Ac__clim, Ac__std = get_climatology_coast('ALK' , SUBLIST, LayerList)
#try:
#    DIC_clim, DIC_std = get_climatology_coast('DIC', SUBLIST, TheMask, LayerList)
#except:
#    print "no adr1 for DIC"
#    DIC_clim, DIC_std = get_climatology_coast('DIC', SUBLIST2, TheMask, LayerList)
#pCO2clim, pCO2std = get_climatology_coast('pCO2',SUBLIST, TheMask, LayerList)
#PH__clim, PH__std = get_climatology_coast('pH' , SUBLIST, TheMask, LayerList)

Ac__clim, Ac__std = get_climatology_coast('ALK' , SUBLIST, LayerList, TheMask)
DIC_clim, DIC_std = get_climatology_coast('DIC', SUBLIST, LayerList, TheMask)
pCO2clim, pCO2std = get_climatology_coast('pCO2',SUBLIST, LayerList, TheMask)
PH__clim, PH__std = get_climatology_coast('pH' , SUBLIST, LayerList, TheMask)


VARLIST=['pCO2','DIC','Ac','pH']
VARLIST=['pCO2','DIC','ALK','pH']
var_dype = [(var,np.float32) for var in VARLIST]
timeseries_DICT={}
for var in VARLIST:
    filename=MODDIR + var + ".pkl"
    TIMESERIES_complete,TL_complete=read_pickle_file(filename)
    ind,ww=TL_complete.select(Req)
    TIMESERIES=TIMESERIES_complete[ind,:]
    timeseries_DICT[var]=TIMESERIES


for iSub, sub in enumerate(basV2.P):
    submask = SubMask(sub, TheMask)
    F = figure_generator.figure_generator(submask)
    fig, axes = F.gen_structure_3(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_4.17_annual." + sub.name + ".png"
    
    for iTime , tTime in enumerate(TL.Timelist):
        datetimelist= [TL.Timelist[iTime] ]
        MEAN = np.zeros((jpk,), dtype=var_dype )
        for ivar, var in enumerate(VARLIST):
            TIMESERIESvar=timeseries_DICT[var]
            mean_profile = TIMESERIESvar[iTime,iSub,icoast,:,istat]
            mean_profile[mean_profile==0]=np.nan
            MEAN[var] = mean_profile
        
        label = 'none'
        color = '0.4'
        figure_generator.profile_plotter(z,MEAN['pCO2'],color, axes[0], None,   label)
        figure_generator.profile_plotter(z,MEAN['DIC' ],color, axes[1], axes[5],label)
#        figure_generator.profile_plotter(z,MEAN['Ac'  ],color, axes[2], axes[6],label)
        figure_generator.profile_plotter(z,MEAN['ALK' ],color, axes[2], axes[6],label)
        figure_generator.profile_plotter(z,MEAN['pH'  ],color, axes[3], axes[7],label)            
        MEAN = np.zeros((jpk,), dtype=var_dype )

    for ivar, var in enumerate(VARLIST):
        TIMESERIESvar=timeseries_DICT[var]
        mean_profile = TIMESERIESvar[:,iSub,icoast,:,istat].mean(axis=0)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile    
    figure_generator.profile_plotter(z,MEAN['pCO2'],'k', axes[0], None,   label)
    figure_generator.profile_plotter(z,MEAN['DIC' ],'k', axes[1], axes[5],label)
#    figure_generator.profile_plotter(z,MEAN['Ac'  ],'k', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['ALK'  ],'k', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['pH'  ],'k', axes[3], axes[7],label) 
    
    Ac__clim_mean=Ac__clim[iSub,:]
    Ac__clim_std =Ac__std[ iSub,:]
    Dic_clim_mean=DIC_clim[iSub,:]
    Dic_clim_std =DIC_std[ iSub,:]
    pCO2clim_mean=pCO2clim[iSub,:]
    pCO2clim_std =pCO2std[ iSub,:]
    PH__clim_mean=PH__clim[iSub,:]
    PH__clim_std =PH__std[ iSub,:]
    
    figure_generator.clim_profile_plotter(z_clim,pCO2clim_mean,pCO2clim_std, axes[0], axes[4])
    figure_generator.clim_profile_plotter(z_clim,Dic_clim_mean,Dic_clim_std, axes[1], axes[5])
    figure_generator.clim_profile_plotter(z_clim,Ac__clim_mean,Ac__clim_std, axes[2], axes[6])
    figure_generator.clim_profile_plotter(z_clim,PH__clim_mean,PH__clim_std, axes[3], axes[7])

    fig.savefig(outfile)
    print (outfile)
    pl.close(fig)

