import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    Produces 3 png files, containing timeseries for some statistics, for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
from commons.mask import Mask
from commons.submask import SubMask

from instruments import lovbio_float as bio_float
from instruments.matchup_manager import Matchup_Manager
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from commons.utils import addsep
from profiler import ALL_PROFILES,TL,BASEDIR
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
from basins.V2 import NRT3 as OGS

def fig_setup(S,subbasin_name):
# ,Lon,Lat):
    from layer_integral import coastline

    fig = plt.figure()
#    ax0 = plt.subplot2grid((4, 3), (0, 0), colspan=2)
#    ax1 = plt.subplot2grid((4, 3), (0, 2))
    ax1 = plt.subplot2grid((4, 3), (0, 0), colspan=3)
    ax2 = plt.subplot2grid((4, 3), (1, 0), colspan=3)
    ax3 = plt.subplot2grid((4, 3), (2, 0), colspan=3)
    ax4 = plt.subplot2grid((4, 3), (3, 0), colspan=3)
    axs = [ax1, ax2, ax3, ax4]
    for ax in [ax2, ax3, ax4]:
        ax.xaxis.grid(True)

    fig.set_size_inches(10,15)
    fig.set_dpi(150)
    c_lon,c_lat=coastline.get()

#    TheMask= Mask(maskfile)
#    S = SubMask(basV2.lev1, maskobject=TheMask)
    bool2d=S.mask_at_level(0)
    smaskplot = np.ones_like(bool2d,dtype=np.float32)
    smaskplot[~bool2d] = np.nan
    lon_min = S.xlevels.min()
    lon_max = S.xlevels.max()
    lat_min = S.ylevels.min()
    lat_max = S.ylevels.max()
    ax1.imshow(smaskplot, extent=[lon_min, lon_max, lat_max, lat_min])
    ax1.invert_yaxis()

    ax1.plot(c_lon,c_lat,'k')
#    ax1.plot(Lon,Lat,'r.')
    ax1.set_title(subbasin_name , color = 'r', fontsize = 18)
    ax1.set_xlim([-10,36])
    ax1.set_ylabel("LAT",color = 'k', fontsize = 15)
    ax1.set_xlabel("LON",color = 'k', fontsize = 15)

#    extent=4
#    ax1.plot(c_lon,c_lat,'k')
#    ax1.plot(Lon,Lat,'ro')
#    ax1.plot(Lon[0],Lat[0],'bo')
#    ax1.set_xlim([Lon.min() -extent/2, Lon.max() +extent/2])
#    ax1.set_ylim([Lat.min() -extent/2, Lat.max() +extent/2])

    return fig, axs

TheMask= Mask(args.maskfile)
INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)
#S = SubMask(basV2.lev1, maskobject=TheMask)

VARLIST = ['P_l','N3n','O2o']
VARLIST_NAME = ['Chlorophyll','Nitrate','Oxygen']
nVar = len(VARLIST)
METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
#MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
#wmo_list=bio_float.get_wmo_list(MED_PROFILES)
#MonthlyRequestors=M.TL.getMonthlist()

times = [req.time_interval.start_time for req in M.TL.getMonthlist()  ]


izmax = TheMask.getDepthIndex(200) 
 
#A_float = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan
A_float = np.load(INDIR + 'Basin_Statistics_FLOAT.npy')
A_model = np.load(INDIR + 'Basin_Statistics_MODEL.npy')

for ivar, var in enumerate(VARLIST):
    for iSub, SubBasin in enumerate(OGS.basin_list):
	OUTFILE = OUTDIR + var + "_" + SubBasin.name + ".png"
        S = SubMask(SubBasin, maskobject=TheMask)
	fig, axes = fig_setup(S,SubBasin.name)
        if (~np.isnan(A_float[ivar,:,iSub,0]).all() == True) or (~np.isnan(A_model[ivar,:,iSub,0]).all() == True):
	    axes[1].plot(times,  A_float[ivar,:,iSub,0],'r',label='REF INTEG')
            axes[1].plot(times,A_model[ivar,:,iSub,0],'b',label='MOD INTEG')

            axes[1].plot(times,A_float[ivar,:,iSub,5],'--r',label='REF SURF')
            axes[1].plot(times,A_model[ivar,:,iSub,5],'--b',label='MOD SURF')
        if (var == "P_l"):
            axes[1].set_ylabel('Chlorophyll \n $[mg{\  } m^{-3}]$',fontsize=15)
        if (var == "O2o"):
            axes[1].set_ylabel('Oxygen 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
        if (var == "N3n"):
            axes[1].set_ylabel('Nitrate 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
        legend = axes[1].legend(loc='upper left', shadow=True, fontsize=12)

        corr = A_float[ivar,:,iSub,1]
#        ref_corr = np.zeros(len(corr), np.float32 ) * np.nan
        for icr , cr in enumerate(corr):
            if (cr <= 0):
                corr[icr] = np.nan
        axes[1].set_xticklabels([])


	if ((np.isnan(corr)).all() == False ): 
	    axes[2].plot(times,corr,'k')
            axes[2].set_ylabel('CORR',fontsize=15)
            axes[2].set_xticklabels([])
	ax2 = axes[2].twinx()
        numP = A_float[ivar,:,iSub,6]
        ax2.bar(times,numP,width=7, color='g', alpha=0.3, label='n points',align='center')
	ax2.set_ylabel(' # Points')
        ax2.set_ylim([0,numP.max() +2])
	ax2.legend(loc=1)
 
        if (var == "P_l"):
            axes[3].plot(times,A_float[ivar,:,iSub,2],'r',label='DCM REF')
            axes[3].plot(times,A_model[ivar,:,iSub,2],'b',label='DCM MOD')
            axes[3].plot(times,A_float[ivar,:,iSub,3],'--r',label='MLB REF')
            axes[3].plot(times,A_model[ivar,:,iSub,3],'--b',label='MLB MOD')

            axes[3].invert_yaxis()
            axes[3].set_ylabel('DCM/MLB $[m]$',fontsize=15)
            axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))

            xlabels = axes[3].get_xticklabels()
            plt.setp(xlabels, rotation=30)
	    legend = axes[3].legend(loc='lower left', shadow=True, fontsize=12)


        if (var == "N3n"):
            axes[3].plot(times,A_float[ivar,:,iSub,4],'r',label='REF')
            axes[3].plot(times,A_model[ivar,:,iSub,4],'b',label='MOD')
            axes[3].invert_yaxis()
            axes[3].set_ylabel('NITRICL $[m]$',fontsize=15)
#        else:
            axes[3].plot(times,  np.ones_like(times))
            axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = axes[3].get_xticklabels()
            plt.setp(xlabels, rotation=30)
            legend = axes[3].legend(loc='upper left', shadow=True, fontsize=12)

        if (var == "O2o"):
            axes[3].plot(times,  np.ones_like(times))
            axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = axes[3].get_xticklabels()
            plt.setp(xlabels, rotation=30)


        fig.savefig(OUTFILE)

# METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf']

#        import sys
#	sys.exit()
