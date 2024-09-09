import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png files for fig4.11 and fig4.19
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
import matplotlib.pyplot as plt
#from commons.layer import Layer
import basins.V2 as OGS
#from static.climatology import get_climatology
#import figure_generator
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from timeseries.plot import Hovmoeller_matrix
from commons.mask import Mask
from commons.submask import SubMask
#import matplotlib.pyplot as pl
from commons.utils import addsep
from profileruns_Summer import runList,colorList


OUTDIR=addsep(args.outdir)
INDIR=addsep(args.inputdir)

TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
TheMask= Mask(args.maskfile)
jpk,jpj,jpi = TheMask.shape
#z = -TheMask.zlevels

#LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]
#z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])


SUBLIST = OGS.Pred.basin_list

    
#N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
#N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList)
#O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList)


VARLIST=['P_l']#,'N1p','N3n','O2o']
#var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

dirSTAT_PROFILES = '/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/'

MEAN = {}
for run in runList:
    print('-------')
    print(run)
    MEAN[run] = np.zeros((72,len(SUBLIST),nVar))
    MODDIR = INDIR + '/RUN_' + run + dirSTAT_PROFILES
    TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")    
    for iSub, sub in enumerate(OGS.Pred):
        print(sub.name)
        submask = SubMask(sub,maskobject=TheMask)
        for ivar,var in enumerate(VARLIST):
            #F = figure_generator.figure_generator(submask)
            #fig, axes = F.gen_structure_1(IDrun,'annual',sub.name)
            #outfile = OUTDIR + "profiles." + sub.name + ".png"  
            Mean_profiles,_,_ = Hovmoeller_matrix(TL.Timelist,TL.filelist, \
                                    var, iSub,coast=1, stat=0, \
                                    depths=np.arange(jpk))
            Mean_profiles[Mean_profiles==0] = np.nan
            MEAN[run][:,iSub,ivar] = np.nanmean(Mean_profiles,1)
                

#plot
print('---------')
print('---------')
print('Plot')
plt.close('all')
fig1,axs = plt.subplots(4,4,num=1,facecolor='w',edgecolor='k', \
                sharex=True,sharey=True) #,figsize=[5,5],dpi=250)

meandiff = np.zeros((len(OGS.Pred.basin_list),len(runList)))
meanreldiff = np.zeros((len(OGS.Pred.basin_list),len(runList)))
for isub,sub in enumerate(OGS.Pred):
    iplot = isub/4
    jplot = isub-iplot*4
    ax1 = axs[iplot,jplot]
    plt.sca(ax1)
    for irun,run in enumerate(runList):
        meandiff[isub,irun] = np.nanmean(np.abs(MEAN['CR'][:,isub,0]-MEAN[run][:,isub,0]))
        meanreldiff[isub,irun] = meandiff[isub,irun]/np.nanmean(MEAN['CR'][:,isub,0])
        plt.plot(MEAN[run][:,isub,0],TheMask.zlevels,'-',color=colorList[run], \
                 label=run)
    plt.ylim(200,0)
    plt.title(sub.name,fontsize=10)
    if (isub==0):
        plt.legend(loc='best',fontsize=8)
    ax1.grid()
    ax1.tick_params(axis='both',labelsize=8)
    ax1.set_xticks([0.,0.1,0.2,0.3])

    fig1.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel('[mg chl/m^3]',fontsize=10)
    plt.ylabel('Depth [m]',fontsize=10)


plt.show(block=False)

outfile = OUTDIR + 'profiles_sub16.png'
plt.savefig(outfile)

for irun,run in enumerate(runList):
    print(irun)
    print(run)
    print(np.nanmean(meandiff[:,irun]))
    print(np.nanmean(meanreldiff[:,irun])*100)
