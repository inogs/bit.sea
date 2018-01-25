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
from instruments import lovbio_float as bio_float
from instruments.matchup_manager import Matchup_Manager
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from commons.utils import addsep
from profiler_floatsat import ALL_PROFILES,TL,RUN
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
import basins.V2 as OGS
from datetime import datetime

def fig_setup(wmo,Lon,Lat):
    from layer_integral import coastline

    fig = plt.figure()
    ax0 = plt.subplot2grid((4, 3), (0, 0), colspan=2)
    ax1 = plt.subplot2grid((4, 3), (0, 2))
    ax2 = plt.subplot2grid((4, 3), (1, 0), colspan=3)
    ax3 = plt.subplot2grid((4, 3), (2, 0), colspan=3)
    ax4 = plt.subplot2grid((4, 3), (3, 0), colspan=3)
    axs = [ax0, ax1, ax2, ax3, ax4]
    for ax in [ax2, ax3, ax4]:
        ax.xaxis.grid(True)
        #ax.set_xlim([datetime(2015,1,1),datetime(2017,1,1)])
        ax.set_xlim(TL.Timelist[0],TL.Timelist[-1])

    fig.set_size_inches(10,15)
    fig.set_dpi(150)
    c_lon,c_lat=coastline.get()

#    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo_list[j])
    ax0.plot(c_lon,c_lat,'k')
    ax0.plot(Lon,Lat,'r.')
    ax0.plot(Lon[0],Lat[0],'b.')
    ax0.set_title("TRAJECTORY of FLOAT " + wmo , color = 'r', fontsize = 18)
#    ind_max_sup=plotmat[0,:].argmax()
    
#    print Lon[ind_max_sup],Lat[ind_max_sup]
#    ax0.plot(Lon[ind_max_sup],Lat[ind_max_sup],'g.')
#    ax0.plot(Lon[0],Lat[0],'bx')
    ax0.set_xlim([-10,36])
    #ax0.set_ylabel("LAT",color = 'k', fontsize = 15)
    #ax0.set_xlabel("LON",color = 'k', fontsize = 15)

    extent=4
    ax1.plot(c_lon,c_lat,'k')
    ax1.plot(Lon,Lat,'ro')
    ax1.plot(Lon[0],Lat[0],'bo')
    ax1.set_xlim([Lon.min() -extent/2, Lon.max() +extent/2])
    ax1.set_ylim([Lat.min() -extent/2, Lat.max() +extent/2])

    return fig, axs


TheMask=Mask(args.maskfile)
INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)

VARLIST = ['P_l']#,'N3n','O2o']
VARLIST_NAME = ['Chlorophyll']#,'Nitrate','Oxygen']
nVar = len(VARLIST)
METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1']
nStat = len(METRICS)

#M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)

izmax = TheMask.getDepthIndex(200) # Max Index for depth 200m

for wmo in wmo_list:
#for wmo in ['6901512']:
    INPUT_FILE = INDIR + wmo + ".nc"
    print INPUT_FILE
    A = ncreader(INPUT_FILE)
    wmo_track_list = bio_float.filter_by_wmo(ALL_PROFILES,wmo)
    nP = len(wmo_track_list)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)
    for ip, p in enumerate(wmo_track_list):
        Lon[ip] = p.lon
        Lat[ip] = p.lat
    times = [p.time for p in wmo_track_list]
    for ivar, var in enumerate(VARLIST):
        OUTFILE = OUTDIR + var + "_" + wmo + ".png"
        print OUTFILE
        fig, axes = fig_setup(wmo,Lon,Lat)
        fig.suptitle(VARLIST_NAME[ivar],fontsize=36,color='b')
        model, ref =A.plotdata(var,'Int_0-200')
        surf_model, surf_ref = A.plotdata(var,'SurfVal')
        if (~np.isnan(model).all() == True) or (~np.isnan(ref).all() == True): 
            axes[2].plot(times,  ref,'r',label='REF INTEG')
            axes[2].plot(times,model,'b',label='MOD INTEG')
            axes[2].plot(times,  surf_ref,'--r',label='REF SURF')
            axes[2].plot(times,  surf_model,'--b',label='MOD SURF')
            if (var == "P_l"):
                axes[2].set_ylabel('Chlorophyll \n $[mg{\  } m^{-3}]$',fontsize=15)
                axes[4].set_ylim(40,200)
            if (var == "O2o"):
                axes[2].set_ylabel('Oxygen 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
            if (var == "N3n"):
                axes[2].set_ylabel('Nitrate 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
                #axes[2].set_ylabel('INTEG 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
                axes[2].set_ylim(0,8)
                axes[4].set_ylim(0,200)
            legend = axes[2].legend(loc='upper left', shadow=True, fontsize=12)
            model_corr, ref_corr =A.plotdata(var,'Corr')
            times_r = times
            for icr , cr in enumerate(ref_corr):
                if (cr <= 0):
                    ref_corr[icr] = np.nan
                    #times_r.remove(times[icr])
            
            axes[3].plot(times,ref_corr,'k')
            #axes[3].plot(times_r,ref_corr[ref_corr>0],'k')
            axes[3].set_ylabel('CORR',fontsize=15)
            axes[2].set_xticklabels([])
            axes[2].set_title(RUN,fontsize=18)
            axes[3].set_ylim(0,1)

        if (var == "P_l"): 
            model_dcm, ref_dcm =A.plotdata(var,'DCM')
            model_mld, ref_mld =A.plotdata(var,'z_01')
            #if (model_mld > 150):
            #model_mld = np.nan
            
            if (~np.isnan(model_dcm).all() == True) or (~np.isnan(ref_dcm).all() == True):
                axes[3].set_xticklabels([])
                axes[4].plot(times,  ref_dcm,'r',label='DCM REF')
                axes[4].plot(times,model_dcm,'b',label='DCM MOD')
                axes[4].plot(times, ref_mld,'--r',label='MLB REF')
                axes[4].plot(times,model_mld,'--b',label='MLB MOD')
                #axes[4].plot(times,np.ones_like(times)* np.nanmean(ref_dcm),'r',linewidth=3) 
                #axes[4].plot(times,np.ones_like(times)* np.nanmean(model_dcm),'b',linewidth=3) #marker='.')
                axes[4].invert_yaxis()
                axes[4].set_ylabel('DCM/MLB $[m]$',fontsize=15)
                #axes[4].xaxis_date()
                #axes[4].xaxis.set_major_locator(mdates.MonthLocator())
                axes[4].xaxis.set_major_formatter(mdates.DateFormatter("%d%b%y"))
                xlabels = axes[4].get_xticklabels()
                plt.setp(xlabels, rotation=30)
            else:
                axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%d%b%y"))
                xlabels = axes[3].get_xticklabels()
                plt.setp(xlabels, rotation=30)

    #	    if (~np.isnan(model_mld).all() == True) or (~np.isnan(ref_mld).all() == True):
    #	        axes_4b = axes[4].twinx()
    #	        axes_4b.plot(times, ref_mld,'--r',label='REF')
    #                axes_4b.plot(times,model_mld,'--b',label='MOD')
    #	        axes_4b.invert_yaxis()
    #	        axes_4b.set_ylabel('MLD $[m]$ - -',fontsize=15)
            legend = axes[4].legend(loc='lower left', shadow=True, fontsize=12)

        if (var == "N3n"):
            model_nit, ref_nit =A.plotdata(var,'Nit_1')
            if (~np.isnan(ref_nit).all() == True) or (~np.isnan(model_nit).all() == True):
            	axes[4].plot(times,  ref_nit,'r',label='REF')
            	axes[4].plot(times,model_nit,'b',label='MOD')
            	axes[4].invert_yaxis()
                axes[4].set_ylabel('NITRICL $[m]$',fontsize=15)
            else: 
		        axes[4].plot(times,  np.ones_like(times))
            axes[4].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = axes[4].get_xticklabels()
            plt.setp(xlabels, rotation=30)
            legend = axes[4].legend(loc='upper left', shadow=True, fontsize=12)
            
        if (var == "O2o"):
            axes[4].plot(times,  np.ones_like(times))
            axes[4].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = axes[4].get_xticklabels()
            plt.setp(xlabels, rotation=30)

        fig.savefig(OUTFILE)
        plt.close(fig)
#    import sys
#    sys.exit()
