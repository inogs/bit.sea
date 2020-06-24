import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png files based on fig4.11
    to compare floats and climatology
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--inwoa', '-w',
                                type = str,
                                default = None,
                                required = True,
                                help = "Directory with files of WOA")
    # parser.add_argument(   '--maskfile', '-m',
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
from instruments import superfloat as bio_float
from static.climatology import get_climatology
from basins.basin import ComposedBasin
from commons.layer import Layer
from commons.mask import Mask
from commons.submask import SubMask
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
from instruments.var_conversions import FLOATVARS
from timeseries.plot import Hovmoeller_matrix
from timeseries.plot import read_pickle_file, read_basic_info

print "start"
IDrun='floatcfr'
OUTDIR=addsep(args.outdir)
DIRWOA=addsep(args.inwoa)
# MODDIR=addsep(args.inputdir)

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
maskfile8="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
Mask8 = Mask(maskfile8)
jpk8,jpj8,jpi8 = Mask8.shape
# TheMask= Mask(args.maskfile, loadtmask=False)
# jpk,jpj,jpi = TheMask.shape
# z = -TheMask.zlevels

PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

# TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")
Profilelist_1=bio_float.FloatSelector(None,TI,basV2.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

SUBLIST = basV2.P.basin_list
nSub = len(SUBLIST)

print 'N3nclim'
N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
print 'N1pclim'
N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList)
print 'O2oclim'
O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList)

# Reading WOA
print '- Reading WOA -'
WOAvar = ['N3n','O2o']
WOAvar = ['N3n','O2o']
DICTwoavar = {
    'N3n': 'n',
    'O2o': 'o',
}

convO2o=1/24.4665E-3 # from BFM 40.8722 con T=25 

woa3D = {}
firstreading = True
for var in WOAvar:
    filewoa = DIRWOA + 'woa13_all_' + DICTwoavar[var] + '00_01.nc'
    print ' ... ' + filewoa
    dataset = netCDF4.Dataset(filewoa)
    if firstreading:
        lonwoa = dataset.variables['lon'][:].data
        latwoa = dataset.variables['lat'][:].data
        levwoa = dataset.variables['depth'][:].data
        levwoa1000 = levwoa[levwoa<=1000]
        nwoa1000 = len(levwoa1000)
        firstreading = False
    woa3D[var] = {}
    for vartype in ['mn','an']:
        varname = DICTwoavar[var] + '_' + vartype
        woa3D[var][vartype] = dataset.variables[varname][0].data
        masknan = woa3D[var][vartype]>10.e+30
        woa3D[var][vartype][masknan] = np.nan
        if var is 'O2o':
            woa3D[var][vartype] = woa3D[var][vartype]*convO2o
    dataset.close()


# VARLIST=['P_l','N1p','N3n','O2o']
VARLIST=['N3n','O2o']
VARLIST=['P_l','N3n','O2o']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

# Adj = {
#     'P_l': True,
#     'N3n': True,
#     'O2o': False,
# }

NewPres_5m=np.linspace(0,1000,201)
nZ = NewPres_5m.shape[0]
zFloat = -NewPres_5m

DICTiSub = {}
for iSub,sub in enumerate(basV2.P):
    DICTiSub[sub.name] = iSub

PresVar = {}
for var_mod in VARLIST:
    PresVar[var_mod] = []
# for wmo in [wmo_list[0]]:
for wmo in wmo_list:
    iip = 0
    print wmo
    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    np_tot = len(list_float_track)
        
    timeseries_DICT = {}
    for var_mod in VARLIST:
        timeseries_DICT[var_mod] = np.zeros((np_tot,nSub,nZ))
        timeseries_DICT[var_mod][:,:,:] = np.nan
  
    sublistall = []
    LIST_indwoa = []
    ipini = -1
    jpini = -1

    for p in list_float_track:
        # for indSub, sub in enumerate(basV2.Pred):
        ip_woa = np.argwhere(p.lon>lonwoa)[-1][0]
        jp_woa = np.argwhere(p.lat>latwoa)[-1][0]
        if (not(ip_woa==ipini)) or (not(jp_woa==jpini)):
            LIST_indwoa.append((jp_woa,ip_woa))
            ipini = ip_woa
            jpini = jp_woa
        for indSub,sub in enumerate(basV2.Pred):
            if sub.is_inside(p.lon,p.lat):
                sublistall.append(sub.name)
                iSubp = indSub

        for var_mod in VARLIST:
            var = FLOATVARS[var_mod]
            # adj=Adj[var_mod]
            try:
                Pres,Prof,Qc=p.read(var)
            except:
                continue
            if Pres.size:
                PresVar[var_mod].append(Pres)
            ii = Pres<=1000 # Deve avere almeno 5 records:
            if len(Prof[ii])>5:
                NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
                timeseries_DICT[var_mod][iip,iSubp,:] = NewProf_5m
                #add all profiles to Med
                # timeseries_DICT[var_mod][iip,-1,:] = NewProf_5m 
        iip += 1
    
    condN3n = np.all(np.isnan(timeseries_DICT['N3n']))
    condO2o = np.all(np.isnan(timeseries_DICT['O2o']))
    if condN3n and condO2o:
        print ' ... Any O2o or N3n profile for ' + wmo
        continue


    sublist = list(np.unique(sublistall))
    subfloatlist = []
    for iSub,sub in enumerate(basV2.P):
        if sub.name in sublist:
            subfloatlist.append(sub)
    
    for iSubf,subf in enumerate(subfloatlist):
        print '.... ' + subf.name
    # subfloat = ComposedBasin(wmo + ' '.join(sublist), \
    #     [sub for sub in subfloatlist],'subfloat for ' + wmo)

# for iSub, sub in enumerate(basV2.P):
    # submask = SubMask(subfloat,maskobject=Mask8)
        submask = SubMask(subf,maskobject=Mask8)
        F = figure_generator.figure_generator(submask)
        fig, axes = F.gen_structure_1(IDrun,'annual',subf.name)
        outfile = OUTDIR + "Fig_float_clim." + wmo + '_' + subf.name + "_nwoa.png"
        
        for iTime in range(np_tot):
            MEAN = np.zeros((nZ,), dtype=var_dype )
            for ivar, var in enumerate(VARLIST):
                TIMESERIES=timeseries_DICT[var]
                mean_profile = TIMESERIES[iTime,DICTiSub[subf.name],:]
                mean_profile[mean_profile==0] = np.nan
                MEAN[var] = mean_profile
            
            label = 'none'
            color = '0.4'
            figure_generator.profile_plotter(zFloat,MEAN['P_l'],color, axes[0], None,   label)
            figure_generator.profile_plotter(zFloat,MEAN['N3n'],color, axes[1], axes[5],label)
            # figure_generator.profile_plotter(zFloat,MEAN['N1p'],color, axes[2], axes[6],label)
            figure_generator.profile_plotter(zFloat,MEAN['O2o'],color, axes[3], axes[7],label)            
            MEAN = np.zeros((nZ,), dtype=var_dype )

        for var in VARLIST:
            TIMESERIES = timeseries_DICT[var]
            mean_profile = np.nanmean(TIMESERIES[:,DICTiSub[subf.name],:],0)
            mean_profile[mean_profile==0]=np.nan
            MEAN[var] = mean_profile
        figure_generator.profile_plotter(zFloat,MEAN['P_l'],'k', axes[0], None,   label)
        figure_generator.profile_plotter(zFloat,MEAN['N3n'],'k', axes[1], axes[5],label)
        # figure_generator.profile_plotter(zFloat,MEAN['N1p'],'k', axes[2], axes[6],label)
        figure_generator.profile_plotter(zFloat,MEAN['O2o'],'k', axes[3], axes[7],label) 

        # for indxswoa in LIST_indwoa:
        #     N3n_woa_mean = woa3D['N3n']['mn'][:nwoa1000,indxswoa[0],indxswoa[1]]
        #     N3n_woa_obja = woa3D['N3n']['an'][:nwoa1000,indxswoa[0],indxswoa[1]]
        #     O2o_woa_mean = woa3D['O2o']['mn'][:nwoa1000,indxswoa[0],indxswoa[1]]
        #     O2o_woa_obja = woa3D['O2o']['an'][:nwoa1000,indxswoa[0],indxswoa[1]]
        
        #     label = 'none' 
        #     figure_generator.profile_plotter(-levwoa1000,N3n_woa_mean,'g', axes[1], axes[5],label)
        #     figure_generator.profile_plotter(-levwoa1000,N3n_woa_obja,'y', axes[1], axes[5],label)
        #     figure_generator.profile_plotter(-levwoa1000,O2o_woa_mean,'g', axes[3], axes[7],label)
        #     figure_generator.profile_plotter(-levwoa1000,O2o_woa_obja,'y', axes[3], axes[7],label)


        N3n_clim_mean=N3n_clim[DICTiSub[subf.name],:]
        N3n_clim_std =N3n_std[DICTiSub[subf.name],:]
        # N1p_clim_mean=N1p_clim[iSub,:]
        # N1p_clim_std =N1p_std[ iSub,:]
        O2o_clim_mean=O2o_clim[DICTiSub[subf.name],:]
        O2o_clim_std =O2o_std[DICTiSub[subf.name],:]
    
        figure_generator.clim_profile_plotter(z_clim,N3n_clim_mean,N3n_clim_std, axes[1], axes[5])
        # figure_generator.clim_profile_plotter(z_clim,N1p_clim_mean,N1p_clim_std, axes[2], axes[6])
        figure_generator.clim_profile_plotter(z_clim,O2o_clim_mean,O2o_clim_std, axes[3], axes[7])

        fig.savefig(outfile)
        print outfile
        pl.close(fig)

