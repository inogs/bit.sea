import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png files of comparison between model and climatology
    in open sea areas (depth higher than 200 m).
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

import matplotlib
matplotlib.use('Agg')

from bitsea.commons.layer import Layer
import bitsea.basins.V2 as basV2
import figure_generator
import figure_generator_extended as fg2
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.timeseries.plot import Hovmoeller_matrix
from bitsea.timeseries.plot import read_pickle_file, read_basic_info
from bitsea.commons import timerequestors
import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
import matplotlib.pyplot as pl
from bitsea.commons.utils import addsep
from bitsea.commons import genUserDateList as DL
import xarray as xr
import os

IDrun='MedBGCins'
OUTDIR=addsep(args.outdir)
MODDIR=addsep(args.inputdir)
CLIMDIR = addsep(args.climdir) + "Clim_Annual_Ref/"

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
Req=timerequestors.Generic_req(TI)

TheMask = Mask.from_file(args.maskfile)
jpk,jpj,jpi = TheMask.shape
z = -TheMask.zlevels

PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")
#weekly=DL.getTimeList(TI.start_time, TI.end_time, days=7)
#TL=TimeList(weekly)

SUBLIST = basV2.Pred.basin_list[:]
SUBLIST.remove(SUBLIST[-1])
VARLIST=['P_l','N1p','N3n','O2o','N4n','N5s']
VARLIST_all=['N1p','N3n','O2o','N4n','N5s','pCO2','DIC','ALK','pH']
# CLIM section for all variables
mean_clim_dict = {}
std_clim_dict = {}
for var in VARLIST_all:
    filename = os.path.join(CLIMDIR, f"{var}_clim_plot.nc")
    with xr.open_dataset(filename) as ds:
        mean_clim_dict[var] = (ds['mean'].dims, ds['mean'].values)
        std_clim_dict[var] = (ds['std'].dims, ds['std'].values)


ds_clim = xr.Dataset(
    {
        f"{var}_clim_mean": xr.DataArray(data, dims=dims)
        for var, (dims, data) in mean_clim_dict.items()
    } | {
        f"{var}_clim_std": xr.DataArray(data, dims=dims)
        for var, (dims, data) in std_clim_dict.items()
    }
)


var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

# reading of input files is limited here ----------
#_, COASTLIST, STAT_LIST = read_basic_info('/g100_scratch/userexternal/vdibiagi/EMODnet_2022/forValidation/ave.20190110-12:00:00.stat_profiles.nc')
_, COASTLIST, STAT_LIST = read_basic_info(TL.filelist[0])
icoast = COASTLIST.index('open_sea')
istat =  STAT_LIST.index('Mean')

timeseries_DICT={}
for var in VARLIST:
    filename=MODDIR + var + ".pkl"
    TIMESERIES_complete,TL_complete=read_pickle_file(filename)
    ind,ww=TL_complete.select(Req)
    TIMESERIES=TIMESERIES_complete[ind,:]
    timeseries_DICT[var]=TIMESERIES
#-------------------------------------------------
var_ordered_nut_clim = ['N3n', 'N1p', 'O2o', 'N5s', 'N4n']
for iSub, sub in enumerate(SUBLIST):
    submask = SubMask(sub, TheMask)
    F = fg2.figure_generator(submask)
    fig, axes = F.gen_structure_1(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_Appendix_nut_open_annual." + sub.name + ".png"

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

    for i, var in enumerate(var_ordered_nut_clim, start=1):#index = 0 is P_l, not plotted as clim
        fg2.clim_profile_plotter(z_clim,ds_clim[f"{var}_clim_mean"][iSub,:],ds_clim[f"{var}_clim_std"][iSub,:],axes[i],axes[i + 6])

    fig.savefig(outfile)
    print (outfile,flush=True)
    pl.close(fig)


###### CARBONATIC VARIABLES


VARLIST=['pCO2','DIC','ALK','pH']
var_dype = [(var,np.float32) for var in VARLIST]
timeseries_DICT={}
for var in VARLIST:
    filename=MODDIR + var + ".pkl"
    TIMESERIES_complete,TL_complete=read_pickle_file(filename)
    ind,ww=TL_complete.select(Req)
    TIMESERIES=TIMESERIES_complete[ind,:]
    timeseries_DICT[var]=TIMESERIES

for iSub, sub in enumerate(SUBLIST):
    submask = SubMask(sub, TheMask)
    F = figure_generator.figure_generator(submask)
    fig, axes = F.gen_structure_3(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_Appendix_carb_open_annual." + sub.name + ".png"

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
    figure_generator.profile_plotter(z,MEAN['ALK'  ],'k', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['pH'  ],'k', axes[3], axes[7],label)
    for i, var in enumerate(VARLIST):
        axis_main = axes[i]
        if i == 0:  # pCO2
            axis_secondary = None
        else:
            axis_secondary = axes[i + 4]
        figure_generator.clim_profile_plotter(z_clim,ds_clim[f"{var}_clim_mean"][iSub,:],ds_clim[f"{var}_clim_std"][iSub,:], axis_main, axis_secondary)



    fig.savefig(outfile)
    print (outfile)
    pl.close(fig)

