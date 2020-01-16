import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png files for fig4.11 and fig4.19
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputdir1','-i1',
                                type = str,
                                required = True,
                                help = '')
    parser.add_argument(   '--inputdir2','-i2',
                                type = str,
                                required = True,
                                help = '')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--maskfile1', '-m1',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--maskfile2', '-m2',
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



from commons.layer import Layer
import basins.V2 as basV2
from static.climatology import get_climatology
from validation.deliverables import figure_generator
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from timeseries.plot import Hovmoeller_matrix
import numpy as np
from commons.mask import Mask
from commons.submask import SubMask
import matplotlib.pyplot as pl
from commons.utils import addsep
IDrun='eas_11'
OUTDIR=addsep(args.outdir)
MODDIR1=addsep(args.inputdir1)
MODDIR2=addsep(args.inputdir2)


TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
maskfile8="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
Mask8 = Mask(maskfile8)
TheMask= Mask(args.maskfile1, loadtmask=False)
jpk_V3,jpj,jpi = TheMask.shape


TheMask_V2= Mask(args.maskfile2, loadtmask=False)
jpk_V2,_,_ = TheMask_V2.shape



PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

TL  = TimeList.fromfilenames(TI, MODDIR1, "ave*nc")
TL2 = TimeList.fromfilenames(TI, MODDIR2, "ave*nc")

SUBLIST = basV2.P.basin_list

    
N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList)
O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList)


VARLIST=['P_l','N1p','N3n','O2o']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)



for iSub, sub in enumerate(basV2.P):
    submask = SubMask(sub,maskobject=Mask8)
    F = figure_generator.figure_generator(submask)
    fig, axes = F.gen_structure_1(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_4.11_annual." + sub.name + ".png"
    
    label="V3"
    jpk = jpk_V3
    MEAN = np.zeros((jpk,), dtype=var_dype )
    z = -TheMask.zlevels
    
    for ivar, var in enumerate(VARLIST):
        Mean_profiles,_,_ = Hovmoeller_matrix(TL.Timelist,TL.filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile  
    figure_generator.profile_plotter(z,MEAN['P_l'],'k', axes[0], None,   label)
    figure_generator.profile_plotter(z,MEAN['N3n'],'k', axes[1], axes[5],label)
    figure_generator.profile_plotter(z,MEAN['N1p'],'k', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['O2o'],'k', axes[3], axes[7],label)

    label="V2"
    jpk = jpk_V2
    MEAN = np.zeros((jpk,), dtype=var_dype )
    z = -TheMask_V2.zlevels

    for ivar, var in enumerate(VARLIST):
        Mean_profiles,_,_ = Hovmoeller_matrix(TL2.Timelist,TL2.filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile
    figure_generator.profile_plotter(z,MEAN['P_l'],'b', axes[0], None,   label)
    figure_generator.profile_plotter(z,MEAN['N3n'],'b', axes[1], axes[5],label)
    figure_generator.profile_plotter(z,MEAN['N1p'],'b', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['O2o'],'b', axes[3], axes[7],label)
        
     
    
    N3n_clim_mean=N3n_clim[iSub,:]
    N3n_clim_std =N3n_std[ iSub,:]
    N1p_clim_mean=N1p_clim[iSub,:]
    N1p_clim_std =N1p_std[ iSub,:]
    O2o_clim_mean=O2o_clim[iSub,:]
    O2o_clim_std =O2o_std[ iSub,:]
    
    figure_generator.clim_profile_plotter(z_clim,N3n_clim_mean,N3n_clim_std, axes[1], axes[5])
    figure_generator.clim_profile_plotter(z_clim,N1p_clim_mean,N1p_clim_std, axes[2], axes[6])
    figure_generator.clim_profile_plotter(z_clim,O2o_clim_mean,O2o_clim_std, axes[3], axes[7])

    figure_generator.add_legend(axes[3])
    fig.savefig(outfile)
    print outfile
    pl.close(fig)



# Figures 4.19
Ac__clim, Ac__std = get_climatology('Ac' , SUBLIST, LayerList)
DIC_clim, DIC_std = get_climatology('DIC', SUBLIST, LayerList)

VARLIST=['pCO2','DIC','Ac','pH']
var_dype = [(var,np.float32) for var in VARLIST]
for iSub, sub in enumerate(basV2.P):
    submask = SubMask(sub,maskobject=Mask8)
    F = figure_generator.figure_generator(submask)
    fig, axes = F.gen_structure_3(IDrun,'annual',sub.name)
    outfile = OUTDIR + "Fig_4.19_annual." + sub.name + ".png"

    label="V3"
    jpk = jpk_V3
    MEAN = np.zeros((jpk,), dtype=var_dype )
    z = -TheMask.zlevels    

    for ivar, var in enumerate(VARLIST):
        Mean_profiles,_,_ = Hovmoeller_matrix(TL.Timelist,TL.filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile    
    figure_generator.profile_plotter(z,MEAN['pCO2'],'k', axes[0], None,   label)
    figure_generator.profile_plotter(z,MEAN['DIC' ],'k', axes[1], axes[5],label)
    figure_generator.profile_plotter(z,MEAN['Ac'  ],'k', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['pH'  ],'k', axes[3], axes[7],label) 

    label="V2"
    jpk = jpk_V2
    MEAN = np.zeros((jpk,), dtype=var_dype )
    z = -TheMask_V2.zlevels
    
    for ivar, var in enumerate(VARLIST):
        Mean_profiles,_,_ = Hovmoeller_matrix(TL2.Timelist,TL2.filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
        mean_profile = Mean_profiles.mean(axis=1)
        mean_profile[mean_profile==0]=np.nan
        MEAN[var] = mean_profile    
    figure_generator.profile_plotter(z,MEAN['pCO2'],'b', axes[0], None,   label)
    figure_generator.profile_plotter(z,MEAN['DIC' ],'b', axes[1], axes[5],label)
    figure_generator.profile_plotter(z,MEAN['Ac'  ],'b', axes[2], axes[6],label)
    figure_generator.profile_plotter(z,MEAN['pH'  ],'b', axes[3], axes[7],label) 

    
    Ac__clim_mean=Ac__clim[iSub,:]
    Ac__clim_std =Ac__std[ iSub,:]
    Dic_clim_mean=DIC_clim[iSub,:]
    Dic_clim_std =DIC_std[ iSub,:]
    
    figure_generator.clim_profile_plotter(z_clim,Dic_clim_mean,Dic_clim_std, axes[1], axes[5])
    figure_generator.clim_profile_plotter(z_clim,Ac__clim_mean,Ac__clim_std, axes[2], axes[6])
    figure_generator.add_legend(axes[3])
    
    fig.savefig(outfile)
    print outfile
    pl.close(fig)

