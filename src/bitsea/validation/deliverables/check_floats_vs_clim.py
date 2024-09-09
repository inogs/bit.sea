import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Exclusion of profiles not consistent with clima at depth
    ''', formatter_class=argparse.RawTextHelpFormatter)

    # parser.add_argument(   '--outdir', '-o',
    #                             type = str,
    #                             default = None,
    #                             required = True,
    #                             help = "")

    # parser.add_argument(   '--inwoa', '-w',
    #                             type = str,
    #                             default = None,
    #                             required = True,
    #                             help = "Directory with files of WOA")
    # # parser.add_argument(   '--maskfile', '-m',
    #                             type = str,
    #                             default = None,
    #                             required = True,
    #                             help = "")
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
import pylab as pl
import basins.V2 as basV2
import figure_generator
import netCDF4
from instruments import lovbio_float as bio_float
from static.climatology import get_climatology
from basins.basin import ComposedBasin
from commons.layer import Layer
from commons.mask import Mask
from commons.submask import SubMask
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
from instruments.var_conversions import LOVFLOATVARS
from timeseries.plot import Hovmoeller_matrix
from timeseries.plot import read_pickle_file, read_basic_info
IDrun='floatcfr'
# OUTDIR=addsep(args.outdir)
# MODDIR=addsep(args.inputdir)

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
# maskfile8="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
# Mask8 = Mask(maskfile8)
# jpk8,jpj8,jpi8 = Mask8.shape
# TheMask= Mask(args.maskfile, loadtmask=False)
# jpk,jpj,jpi = TheMask.shape
# z = -TheMask.zlevels

PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

z_clim = np.array([(l.bottom+l.top)/2  for l in LayerList])

# TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")
Profilelist_1=bio_float.FloatSelector(None,TI,basV2.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

SUBLIST = basV2.P.basin_list
nSub = len(SUBLIST)


N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList)
O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList)



# VARLIST=['P_l','N1p','N3n','O2o']
VARLIST=['N3n']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

Adj = {
    'P_l': True,
    'N3n': True,
    'O2o': False,
}

NewPres_5m=np.linspace(0,40,9)
NewPres_10m=np.linspace(50,200,(200-50)/10+1)
NewPres_25m=np.linspace(225,1000,(1000-225)/25+1)

NewPres = np.concatenate([NewPres_5m,NewPres_10m,NewPres_25m])
maskPres = (NewPres>600) & (NewPres<800)

nZ = NewPres.shape[0]
zFloat = -NewPres

DICTiSub = {}
for iSub,sub in enumerate(basV2.P):
    DICTiSub[sub.name] = iSub

PresVar = {}
for var_mod in VARLIST:
    PresVar[var_mod] = []

LIST = {}
for var_mod in VARLIST:
    LIST[var_mod] = [[] for ii in range(3)]

thresholdRMSD = 1.5

# for wmo in [wmo_list[0]]:
for wmo in wmo_list:
    iip = {}
    iipcheck = {}
    for var_mod in VARLIST:
        iip[var_mod] = 0
        iipcheck[var_mod] = 0
    print wmo
    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    np_tot = len(list_float_track)
    print '  number of profiles for this float ' + np.str(np_tot) 
        
    timeseries_DICT = {}
    for var_mod in VARLIST:
        timeseries_DICT[var_mod] = np.zeros((np_tot,nSub,nZ))
        timeseries_DICT[var_mod][:,:,:] = np.nan
  

    for p in list_float_track:
        # for indSub, sub in enumerate(basV2.Pred):
        for indSub,sub in enumerate(basV2.Pred):
            if sub.is_inside(p.lon,p.lat):
                iSubp = indSub

        for var_mod in VARLIST:
            var = LOVFLOATVARS[var_mod]
            adj=Adj[var_mod]
            Pres,Prof,Qc=p.read(var,read_adjusted=adj)
            if Pres.size:
                PresVar[var_mod].append(Pres)
            ii = Pres<=1000 # Deve avere almeno 5 records:
            if len(Prof[ii])>5:
                iip[var_mod] += 1
                NewProf = np.interp(NewPres,Pres[ii],Prof[ii])
                NewClim = np.interp(NewPres,z_clim,N3n_clim[iSubp,:])
                RMSDmasked = (np.nanmean((NewProf[maskPres]-NewClim[maskPres])**2))**.5
                if RMSDmasked>thresholdRMSD :
                    # print '    ...  excluding a profile for comparison with clim'
                    iipcheck[var_mod] += 1
                # timeseries_DICT[var_mod][iip,iSubp,:] = NewProf
                #add all profiles to Med
                # timeseries_DICT[var_mod][iip,-1,:] = NewProf_5m 


    for var_mod in VARLIST:
        if iip[var_mod]==0:
            print '    Any profile for ' + var_mod
        else:
            LIST[var_mod][0].append(wmo[-3:])
            LIST[var_mod][1].append(iip[var_mod])
            LIST[var_mod][2].append(iip[var_mod]-iipcheck[var_mod])

            percex = np.float(iipcheck[var_mod])/np.float(iip[var_mod])*100.
            print '  Number of profiles for ' + var_mod + \
                ' ' + np.str(iip[var_mod])
            print '  Number of excluded ' + \
                '  ' + np.str(iipcheck[var_mod]) + ' percent ' + \
                np.str(percex) + '  (considered ' + np.str(100.-percex) + ')'


pl.close('all')
fig,axs = pl.subplots(3,1,figsize=[12,8],sharex=True)

pl.sca(axs[0])
pl.title('Number of profiles available for N3n')
pl.bar(LIST['N3n'][0],LIST['N3n'][1])
pl.grid()

pl.sca(axs[1])
pl.title('Number of non-excluded profiles')
pl.bar(LIST['N3n'][0],LIST['N3n'][2])
pl.grid()

pl.sca(axs[2])
pl.title('Percentage of non-excluded profiles')
percnoex = np.array(LIST['N3n'][2]).astype(float) / np.array(LIST['N3n'][1]).astype(float) * 100.
pl.bar(LIST['N3n'][0],percnoex)
pl.grid()

pl.show(block=False)
pl.savefig('excluded_profiles_2015.png')
    
#     condN3n = np.all(np.isnan(timeseries_DICT['N3n']))
#     condO2o = np.all(np.isnan(timeseries_DICT['O2o']))
#     if condN3n and condO2o:
#         print ' ... Any O2o or N3n profile for ' + wmo
#         continue


    
#     for iSubf,subf in enumerate(subfloatlist):
#         print '.... ' + subf.name
#     # subfloat = ComposedBasin(wmo + ' '.join(sublist), \
#     #     [sub for sub in subfloatlist],'subfloat for ' + wmo)

# # for iSub, sub in enumerate(basV2.P):
#     # submask = SubMask(subfloat,maskobject=Mask8)
#         submask = SubMask(subf,maskobject=Mask8)
#         F = figure_generator.figure_generator(submask)
#         fig, axes = F.gen_structure_1(IDrun,'annual',subf.name)
#         outfile = OUTDIR + "Fig_float_clim." + wmo + '_' + subf.name + ".png"
        
#         for iTime in range(np_tot):
#             MEAN = np.zeros((nZ,), dtype=var_dype )
#             for ivar, var in enumerate(VARLIST):
#                 TIMESERIES=timeseries_DICT[var]
#                 mean_profile = TIMESERIES[iTime,DICTiSub[subf.name],:]
#                 mean_profile[mean_profile==0] = np.nan
#                 MEAN[var] = mean_profile
            
#             label = 'none'
#             color = '0.4'
#             figure_generator.profile_plotter(zFloat,MEAN['P_l'],color, axes[0], None,   label)
#             figure_generator.profile_plotter(zFloat,MEAN['N3n'],color, axes[1], axes[5],label)
#             # figure_generator.profile_plotter(zFloat,MEAN['N1p'],color, axes[2], axes[6],label)
#             figure_generator.profile_plotter(zFloat,MEAN['O2o'],color, axes[3], axes[7],label)            
#             MEAN = np.zeros((nZ,), dtype=var_dype )

#         for var in VARLIST:
#             TIMESERIES = timeseries_DICT[var]
#             mean_profile = np.nanmean(TIMESERIES[:,DICTiSub[subf.name],:],0)
#             mean_profile[mean_profile==0]=np.nan
#             MEAN[var] = mean_profile
#         figure_generator.profile_plotter(zFloat,MEAN['P_l'],'k', axes[0], None,   label)
#         figure_generator.profile_plotter(zFloat,MEAN['N3n'],'k', axes[1], axes[5],label)
#         # figure_generator.profile_plotter(zFloat,MEAN['N1p'],'k', axes[2], axes[6],label)
#         figure_generator.profile_plotter(zFloat,MEAN['O2o'],'k', axes[3], axes[7],label) 

#         for indxswoa in LIST_indwoa:
#             N3n_woa_mean = woa3D['N3n']['mn'][:nwoa1000,indxswoa[0],indxswoa[1]]
#             N3n_woa_obja = woa3D['N3n']['an'][:nwoa1000,indxswoa[0],indxswoa[1]]
#             O2o_woa_mean = woa3D['O2o']['mn'][:nwoa1000,indxswoa[0],indxswoa[1]]
#             O2o_woa_obja = woa3D['O2o']['an'][:nwoa1000,indxswoa[0],indxswoa[1]]
        
#             label = 'none' 
#             figure_generator.profile_plotter(-levwoa1000,N3n_woa_mean,'g', axes[1], axes[5],label)
#             figure_generator.profile_plotter(-levwoa1000,N3n_woa_obja,'y', axes[1], axes[5],label)
#             figure_generator.profile_plotter(-levwoa1000,O2o_woa_mean,'g', axes[3], axes[7],label)
#             figure_generator.profile_plotter(-levwoa1000,O2o_woa_obja,'y', axes[3], axes[7],label)


#         N3n_clim_mean=N3n_clim[DICTiSub[subf.name],:]
#         N3n_clim_std =N3n_std[DICTiSub[subf.name],:]
#         # N1p_clim_mean=N1p_clim[iSub,:]
#         # N1p_clim_std =N1p_std[ iSub,:]
#         O2o_clim_mean=O2o_clim[DICTiSub[subf.name],:]
#         O2o_clim_std =O2o_std[DICTiSub[subf.name],:]
    
#         figure_generator.clim_profile_plotter(z_clim,N3n_clim_mean,N3n_clim_std, axes[1], axes[5])
#         # figure_generator.clim_profile_plotter(z_clim,N1p_clim_mean,N1p_clim_std, axes[2], axes[6])
#         figure_generator.clim_profile_plotter(z_clim,O2o_clim_mean,O2o_clim_std, axes[3], axes[7])

#         fig.savefig(outfile)
#         print outfile
#         pl.close(fig)

