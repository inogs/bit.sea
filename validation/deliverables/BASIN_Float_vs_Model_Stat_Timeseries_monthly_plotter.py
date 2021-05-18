import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    It reads in the input directory two files ( model and ref) 
    containing [(nVar, nTime, nSub, nStat)] arrays to generate
    the following metrics:

    CHL-PROF-D-CLASS4-PROF-CORR-BASIN
    NIT-PROF-D-CLASS4-PROF-CORR-BASIN
     DO-PROF-D-CLASS4-PROF-CORR-BASIN 

    FIGURES    - It produces 3 png files, containing timeseries for some statistics, for each basin (7 basins).
    FIG.IV.5   - for Chla
    FIG.IV.13  - for Nitrate

    TABLES     - It produces also a table for the relative MEANS, BIAS and RMSD.
    TABLE IV.4 - for Chla
    TABLE IV.7 - for Nitrate
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
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

from instruments import superfloat as bio_float
from instruments.matchup_manager import Matchup_Manager
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from commons.utils import addsep
from profiler_floats import ALL_PROFILES,TL,BASEDIR
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
from basins.V2 import NRT3 as OGS
from commons.time_interval import TimeInterval
from matchup.statistics import *
from commons.utils import writetable
from datetime import datetime
from datetime import timedelta
from profiler_floats import *

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
#        ax.set_xlim([datetime(2015,1,1),datetime(2017,1,1)])
        ax.set_xlim([datetime(int(DATESTART[0:4]),int(DATESTART[4:6]),int(DATESTART[6:9])),datetime(int(DATE__END[0:4]),int(DATE__END[4:6]),int(DATE__END[6:9]))+timedelta(days=1)])

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

VARLIST      = ['P_l','N3n','O2o','P_c']
VARLIST_NAME = ['Chlorophyll','Nitrate','Oxygen','PhytoC']
nVar         = len(VARLIST)
METRICS      = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf']


METRICS_ALL=['CORR','INTmeanRef','INTmeanMod','INTstdRef','INTstdMod','INT_0-200_BIAS','INT_0-200_RMSD','DCMmeanRef','DCMmeanMod','DCMstdRef','DCMstdMod','DCM_BIAS','DCM_RMSD','WBLmeanRef','WBLmeanMod','WBLstdRef','WBLstdMod','WBL_BIAS','WBL_RMSD','NITmeanRef','NITmeanMod','NITstdRef','NITstdMod','NIT_BIAS','NIT_RMSD','OMZmeanRef','OMZmeanMod','OMZstdRef','OMZstdMod','OMZ_BIAS','OMZ_RMSD','O2oMaxRef','O2oMaxMod','O2oMaxRef_std','O2oMaxMod_std','O2oMax_BIAS','O2oMax_RMSD','N_POINTS']

nm=len(METRICS_ALL)
print nm
METRICS_SHORT= ['CORR','INT_0-200_BIAS','INT_0-200_RMSD','DCM_BIAS','DCM_RMSD','WBL_BIAS','WBL_RMSD','NIT_BIAS','NIT_RMSD','OMZ_BIAS','OMZ_RMSD','O2oMax_BIAS','O2oMax_RMSD','N_POINTS']
nStat        = len(METRICS)
nSub         = len(OGS.basin_list)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
#MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
#wmo_list=bio_float.get_wmo_list(MED_PROFILES)
#MonthlyRequestors=M.TL.getMonthlist()

times = [req.time_interval.start_time for req in M.TL.getMonthlist()  ]
ti_restrict = TimeInterval(DATESTART,DATE__END,"%Y%m%d")
ii = np.zeros((len(times),) , np.bool)
for k,t in enumerate(times) : ii[k] = ti_restrict.contains(t)


izmax = TheMask.getDepthIndex(200) 
 
#A_float = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan
A_float = np.load(INDIR + 'Basin_Statistics_FLOAT.npy')
A_model = np.load(INDIR + 'Basin_Statistics_MODEL.npy')

for ivar, var in enumerate(VARLIST):
# if ( var == "P_c"):
    TABLE_METRICS = np.zeros((nSub,nm),np.float32)*np.nan #before 18 metrics
    TABLE_METRICS_SHORT = np.zeros((nSub,14),np.float32)*np.nan #before 10 metrics
    for iSub, SubBasin in enumerate(OGS.basin_list):
	OUTFILE = OUTDIR + var + "_" + SubBasin.name + ".png"
        S = SubMask(SubBasin, maskobject=TheMask)
	fig, axes = fig_setup(S,SubBasin.name)
        if (~np.isnan(A_float[ivar,:,iSub,0]).all() == True) or (~np.isnan(A_model[ivar,:,iSub,0]).all() == True):
	    Int_Ref = A_float[ivar,:,iSub,0]
	    Int_Mod = A_model[ivar,:,iSub,0]
	    axes[1].plot(times, Int_Ref,'r',label='REF INTEG')
            axes[1].plot(times,Int_Mod,'b',label='MOD INTEG')
            
	    TABLE_METRICS[iSub,1] = np.nanmean(Int_Ref[ii])
	    TABLE_METRICS[iSub,2] = np.nanmean(Int_Mod[ii])
            TABLE_METRICS[iSub,3] = np.nanstd(Int_Ref[ii])
            TABLE_METRICS[iSub,4] = np.nanstd(Int_Mod[ii])

            good = ~np.isnan(Int_Ref[ii]) &  ~np.isnan(Int_Mod[ii])
            m = matchup(Int_Mod[ii][good],   Int_Ref[ii][good])

            TABLE_METRICS[iSub,5] = m.bias()
            TABLE_METRICS[iSub,6] = m.RMSE()

             
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
	    axes[2].set_ylim([0.2,1])
            TABLE_METRICS[iSub,0] = np.nanmean(corr[ii])

	ax2 = axes[2].twinx()
        numP = A_float[ivar,:,iSub,6]
        ax2.bar(times,numP,width=7, color='g', alpha=0.3, label='n points',align='center')
	ax2.set_ylabel(' # Points')
        ax2.set_ylim([0,numP.max() +2])
	ax2.legend(loc=1)
	TABLE_METRICS[iSub,37] = np.mean(numP[ii])
 
        if (var == "P_l"):
	    DCM_Ref = A_float[ivar,:,iSub,2]
	    DCM_Mod = A_model[ivar,:,iSub,2]
	    WBL_Ref = A_float[ivar,:,iSub,3]
	    WBL_Mod = A_model[ivar,:,iSub,3]
            axes[3].plot(times,DCM_Ref,'r',label='DCM REF')
            axes[3].plot(times,DCM_Mod,'b',label='DCM MOD')
	    # FILTER OUT MAY TO NOV INCLUDED:
	    WBL_Ref[4:11] = np.nan
	    WBL_Mod[4:11] = np.nan
            WBL_Ref[16:23] = np.nan
            WBL_Mod[16:23] = np.nan
            axes[3].plot(times,WBL_Ref,'--r',label='WBL REF')
            axes[3].plot(times,WBL_Mod,'--b',label='WBL MOD')

            axes[3].invert_yaxis()
            axes[3].set_ylabel('DCM/WBL $[m]$',fontsize=15)
            axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            axes[1].set_ylim([0,1.0])
	    ax2.set_ylim([0,50])
            axes[3].set_ylim([160,40])

            xlabels = axes[3].get_xticklabels()
            plt.setp(xlabels, rotation=30)
	    legend = axes[3].legend(loc='lower right', shadow=True, fontsize=12)

            TABLE_METRICS[iSub,7] = np.nanmean(DCM_Ref[ii])
            TABLE_METRICS[iSub,8] = np.nanmean(DCM_Mod[ii])
            TABLE_METRICS[iSub,9] = np.nanstd(DCM_Ref[ii])
            TABLE_METRICS[iSub,10] = np.nanstd(DCM_Mod[ii])

            good = ~np.isnan(DCM_Ref[ii]) &  ~np.isnan(DCM_Mod[ii])
            m = matchup(DCM_Mod[ii][good],   DCM_Ref[ii][good])
            TABLE_METRICS[iSub,11] = m.bias()
            TABLE_METRICS[iSub,12] = m.RMSE()

            TABLE_METRICS[iSub,13] = np.nanmean(WBL_Ref[ii])
            TABLE_METRICS[iSub,14] = np.nanmean(WBL_Mod[ii])
            TABLE_METRICS[iSub,15] = np.nanstd(WBL_Ref[ii])
            TABLE_METRICS[iSub,16] = np.nanstd(WBL_Mod[ii])

            good = ~np.isnan(WBL_Ref[ii]) &  ~np.isnan(WBL_Mod[ii])
            m = matchup(WBL_Mod[ii][good],   WBL_Ref[ii][good])
            TABLE_METRICS[iSub,17] = m.bias()
            TABLE_METRICS[iSub,18] = m.RMSE()



        if (var == "N3n"):
	    Nit_Ref = A_float[ivar,:,iSub,4]
	    Nit_Mod = A_model[ivar,:,iSub,4]
            axes[3].plot(times,Nit_Ref,'r',label='REF')
            axes[3].plot(times,Nit_Mod,'b',label='MOD')
            axes[3].invert_yaxis()
            axes[3].set_ylabel('NITRACL $[m]$',fontsize=15)

            TABLE_METRICS[iSub,19] = np.nanmean(Nit_Ref[ii])
            TABLE_METRICS[iSub,20] = np.nanmean(Nit_Mod[ii])
            TABLE_METRICS[iSub,21] = np.nanstd(Nit_Ref[ii])
            TABLE_METRICS[iSub,22] = np.nanstd(Nit_Mod[ii])
            good = ~np.isnan(Nit_Ref[ii]) &  ~np.isnan(Nit_Mod[ii])
            m = matchup(Nit_Mod[ii][good],   Nit_Ref[ii][good])
	    TABLE_METRICS[iSub,23] = m.bias()
            TABLE_METRICS[iSub,24] = m.RMSE()
#        else:
            axes[3].plot(times,  np.ones_like(times))
            axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = axes[3].get_xticklabels()
            plt.setp(xlabels, rotation=30)
            legend = axes[3].legend(loc='lower right', shadow=True, fontsize=12)
            axes[1].set_ylim([0,6])
            ax2.set_ylim([0,45])
            axes[3].set_ylim([200,0])

        if (var == "O2o"):
#            axes[3].plot(times,  np.ones_like(times))
            OMZmeanRef = A_float[ivar,:,iSub,7]
            OMZmeanMod = A_model[ivar,:,iSub,7]
            O2oMaxRef = A_float[ivar,:,iSub,8]
            O2oMaxMod = A_model[ivar,:,iSub,8]
            axes[3].plot(times,O2oMaxRef,'r',label='REF')
            axes[3].plot(times,O2oMaxMod,'b',label='MOD')
            axes[3].invert_yaxis()
            axes[3].set_ylabel('MAX OXY depth $[m]$',fontsize=15)
            axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
	    legend = axes[3].legend(loc='lower right', shadow=True, fontsize=12)
            xlabels = axes[3].get_xticklabels()
            plt.setp(xlabels, rotation=30)

            TABLE_METRICS[iSub,25] = np.nanmean(OMZmeanRef[ii])
            TABLE_METRICS[iSub,26] = np.nanmean(OMZmeanMod[ii])
            TABLE_METRICS[iSub,27] = np.nanstd(OMZmeanRef[ii])
            TABLE_METRICS[iSub,28] = np.nanstd(OMZmeanMod[ii])

            good = ~np.isnan(OMZmeanRef[ii]) &  ~np.isnan(OMZmeanMod[ii])
            m = matchup(OMZmeanMod[ii][good],   OMZmeanRef[ii][good])
            TABLE_METRICS[iSub,29] = m.bias()
            TABLE_METRICS[iSub,30] = m.RMSE()

            TABLE_METRICS[iSub,31] = np.nanmean(O2oMaxRef[ii])
            TABLE_METRICS[iSub,32] = np.nanmean(O2oMaxMod[ii])
            TABLE_METRICS[iSub,33] = np.nanstd(O2oMaxRef[ii])
            TABLE_METRICS[iSub,34] = np.nanstd(O2oMaxMod[ii])

            good = ~np.isnan(O2oMaxRef[ii]) &  ~np.isnan(O2oMaxMod[ii])
            m = matchup(O2oMaxMod[ii][good],   O2oMaxRef[ii][good])
            TABLE_METRICS[iSub,35] = m.bias()
            TABLE_METRICS[iSub,36] = m.RMSE()

        if ( (np.isnan(corr)).all() == True ) : continue
        fig.savefig(OUTFILE)
	plt.close(fig)

    row_names   =[sub.name for sub in OGS.basin_list]

    TABLE_METRICS_SHORT[:,0] = TABLE_METRICS[:,0]
    TABLE_METRICS_SHORT[:,1] = TABLE_METRICS[:,5]
    TABLE_METRICS_SHORT[:,2] = TABLE_METRICS[:,6]
    TABLE_METRICS_SHORT[:,3] = TABLE_METRICS[:,11]
    TABLE_METRICS_SHORT[:,4] = TABLE_METRICS[:,12]
    TABLE_METRICS_SHORT[:,5] = TABLE_METRICS[:,17]
    TABLE_METRICS_SHORT[:,6] = TABLE_METRICS[:,18]
    TABLE_METRICS_SHORT[:,7] = TABLE_METRICS[:,23]
    TABLE_METRICS_SHORT[:,8] = TABLE_METRICS[:,24]
    TABLE_METRICS_SHORT[:,9] = TABLE_METRICS[:,29]

    TABLE_METRICS_SHORT[:,10] = TABLE_METRICS[:,30]
    TABLE_METRICS_SHORT[:,11] = TABLE_METRICS[:,35]
    TABLE_METRICS_SHORT[:,12] = TABLE_METRICS[:,36]
    TABLE_METRICS_SHORT[:,13] = TABLE_METRICS[:,37]


    writetable(OUTDIR + var + '_tab_statistics_ALL.txt',TABLE_METRICS,row_names,METRICS_ALL,fmt="%3.2f\t %3.2f\t %3.2f\t %3.2f\t %3.2f\t %3.2f\t %3.2f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t%.0f\t %.0f\t %.0f\t %.0f\t")

    writetable(OUTDIR + var + '_tab_statistics_SHORT.txt',TABLE_METRICS_SHORT,row_names,METRICS_SHORT,fmt="%3.2f\t %3.2f\t %3.2f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f %.0f\t %.0f\t %.0f\t %.0f\t")

