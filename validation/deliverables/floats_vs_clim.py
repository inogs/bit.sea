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
from instruments import lovbio_float as bio_float
from static.climatology import get_climatology
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
OUTDIR=addsep(args.outdir)
# MODDIR=addsep(args.inputdir)

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
maskfile8="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
Mask8 = Mask(maskfile8)
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


N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList)
O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList)


# VARLIST=['P_l','N1p','N3n','O2o']
VARLIST=['P_l','N3n']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

Adj = {
    'P_l': True,
    'N3n': True,
}

# reading of input files is limited here ----------
# _, COASTLIST, STAT_LIST = read_basic_info(TL.filelist[0])
# icoast = COASTLIST.index('open_sea')
# istat =  STAT_LIST.index('Mean')

# timeseries_DICT={}
# for var in VARLIST:
#     filename=MODDIR + var + ".pkl"
#     TIMESERIES,_=read_pickle_file(filename)
#     timeseries_DICT[var]=TIMESERIES
#-------------------------------------------------

np_tot = 0
for wmo in wmo_list:
# for j in range(0,1):
    list_float_track = bio_float.filter_by_wmo(Profilelist_1,wmo)
    nP = len(list_float_track)
    np_tot = np_tot + nP

NewPres_5m=np.linspace(0,1000,201)
nZ = NewPres_5m.shape[0]
zFloat = -NewPres_5m

timeseries_DICT = {}
for var_mod in VARLIST:
    timeseries_DICT[var_mod] = np.zeros((np_tot,nSub,nZ))
    timeseries_DICT[var_mod][:,:,:] = np.nan

iip = 0
for wmo in wmo_list:
    print np.str(iip) + ' of ' + np.str(np_tot)
    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    print wmo
  
    for p in list_float_track:
        for indSub, sub in enumerate(basV2.Pred):
            if sub.is_inside(p.lon,p.lat):
                iSub = indSub

        for var_mod in VARLIST:
            var = LOVFLOATVARS[var_mod]
            adj=Adj[var_mod]
            Pres,Prof,Qc=p.read(var,read_adjusted=adj)
            ii = Pres<=1000 # Deve avere almeno 5 records:
            if len(Prof[ii])>5:
                NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
                timeseries_DICT[var_mod][iip,iSub,:] = NewProf_5m
                #add all profiles to Med
                timeseries_DICT[var_mod][iip,-1,:] = NewProf_5m 
        iip += 1    

for iSub, sub in enumerate(basV2.P):
    submask = SubMask(sub,maskobject=Mask8)
    F = figure_generator.figure_generator(submask)
    fig, axes = F.gen_structure_1(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_float_clim." + sub.name + ".png"
    
    for iTime in range(np_tot):
        # datetimelist= [TL.Timelist[iTime] ]
        MEAN = np.zeros((nZ,), dtype=var_dype )
        for ivar, var in enumerate(VARLIST):
            TIMESERIES=timeseries_DICT[var]
            #Mean_profiles,_,_ = Hovmoeller_matrix(datetimelist,[filename], var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
            #mean_profile = Mean_profiles[:,0]#.mean(axis=1)
            mean_profile = TIMESERIES[iTime,iSub,:]
            mean_profile[mean_profile==0] = np.nan
            MEAN[var] = mean_profile
        
        label = 'none'
        color = '0.4'
        figure_generator.profile_plotter(zFloat,MEAN['P_l'],color, axes[0], None,   label)
        figure_generator.profile_plotter(zFloat,MEAN['N3n'],color, axes[1], axes[5],label)
        # figure_generator.profile_plotter(zFloat,MEAN['N1p'],color, axes[2], axes[6],label)
        # figure_generator.profile_plotter(zFloat,MEAN['O2o'],color, axes[3], axes[7],label)            
        MEAN = np.zeros((nZ,), dtype=var_dype )

    for var in VARLIST:
        TIMESERIES = timeseries_DICT[var]
        mean_profile = np.nanmean(TIMESERIES[:,iSub,:],0)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile
    figure_generator.profile_plotter(zFloat,MEAN['P_l'],'k', axes[0], None,   label)
    figure_generator.profile_plotter(zFloat,MEAN['N3n'],'k', axes[1], axes[5],label)
    # figure_generator.profile_plotter(zFloat,MEAN['N1p'],'k', axes[2], axes[6],label)
    # figure_generator.profile_plotter(zFloat,MEAN['O2o'],'k', axes[3], axes[7],label) 
    
    N3n_clim_mean=N3n_clim[iSub,:]
    N3n_clim_std =N3n_std[ iSub,:]
    # N1p_clim_mean=N1p_clim[iSub,:]
    # N1p_clim_std =N1p_std[ iSub,:]
    # O2o_clim_mean=O2o_clim[iSub,:]
    # O2o_clim_std =O2o_std[ iSub,:]
    
    figure_generator.clim_profile_plotter(z_clim,N3n_clim_mean,N3n_clim_std, axes[1], axes[5])
    # figure_generator.clim_profile_plotter(z_clim,N1p_clim_mean,N1p_clim_std, axes[2], axes[6])
    # figure_generator.clim_profile_plotter(z_clim,O2o_clim_mean,O2o_clim_std, axes[3], axes[7])

    fig.savefig(outfile)
    print outfile
    pl.close(fig)

