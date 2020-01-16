
import numpy as np
import matplotlib.pyplot as pl
import argparse
#import pickle
import matplotlib.dates as mdates

from basins import OGS
from instruments import lovbio_float as bio_float
from profileruns import runList,colorList
from profiler_floatsat import ALL_PROFILES,TL
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
from commons.utils import writetable

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required =True,
                            default = "DA_FLOAT_SAT/Winter/",
                            help = ''' Input dir of .pkl files produced by ScMYvalidation'''
                            )

    parser.add_argument(   '--outdir', '-O',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )
    # parser.add_argument(   '--open_sea_file', '-o',
    #                         type = str,
    #                         required = True,
    #                         help = 'Input pickle file for open sea')

    return parser.parse_args()

args = argument()

VARLIST = ['P_l']
nVar = len(VARLIST)

MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)
#############

matstats = {}
for var in VARLIST:
    matstats[var] = {}
    for run in runList:
        matstats[var][run] = np.zeros((len(wmo_list),4))

wmo_valid_list = {}
for var in VARLIST:
    wmo_valid_list[var] = []

for iwmo,wmo in enumerate(wmo_list):
# for wmo in ['6901510']:
    print('Float ' + wmo)
    A = {}
    print('Reading input ----')
    for irun,run in enumerate(runList):
        INDIR = args.inputdir + '/RUN_' + run + '/tmp_nc/'
        INPUT_FILE = INDIR + wmo + ".nc"
        A[run] = ncreader(INPUT_FILE)
        wmo_track_list = bio_float.filter_by_wmo(ALL_PROFILES,wmo)
        nP = len(wmo_track_list)
        times = [p.time for p in wmo_track_list]
    
    for var in VARLIST:
        pl.close('all')
        fig,axes = pl.subplots(2,1,sharex=True,dpi=150)
        fig.suptitle(var + ' - ' + wmo,fontsize=12,color='b')
        OUTFILE = args.outdir + var + "_" + wmo + ".png"
        for irun,run in enumerate(runList):
            print('Run ' + run)
            model, ref = A[run].plotdata(var,'Int_0-200')
            surf_model, surf_ref = A[run].plotdata(var,'SurfVal')
            if (~np.isnan(model).all() == True) or (~np.isnan(ref).all() == True):
                matstats[var][run][iwmo,0] = (np.nanmean((model-ref)**2))**0.5
                matstats[var][run][iwmo,1] = np.nanmean(model-ref)
                matstats[var][run][iwmo,2] = (np.nanmean((surf_model-surf_ref)**2))**0.5
                matstats[var][run][iwmo,3] = np.nanmean(surf_model-surf_ref)
                if irun==0:
                    wmo_valid_list[var].append(wmo)
                    axes[0].plot(times,ref,'o-g',label='Float')
                axes[0].plot(times,model,'-',color=colorList[run],label=run)
                axes[0].set_title('Integral 0-200m',fontsize=10)

                if irun==0:
                    axes[1].plot(times,surf_ref,'o-g',label='Float')
                axes[1].plot(times,surf_model,'-',color=colorList[run],label=run)
                axes[1].set_title('Surface',fontsize=10)
                if (var == "P_l"):
                    axes[0].set_ylabel('Chlorophyll \n $[mg{\  } m^{-3}]$',fontsize=10)
                    axes[1].set_ylabel('Chlorophyll \n $[mg{\  } m^{-3}]$',fontsize=10)
                if (var == "O2o"):
                    axes[0].set_ylabel('Oxygen 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=10)
                if (var == "N3n"):
                    axes[0].set_ylabel('Nitrate 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=10)
        

        for iax in range(2):
            axes[iax].grid(True)
            legend0 = axes[iax].legend(loc='best', fontsize=8)
            axes[iax].tick_params(axis='both',labelsize=8)
        axes[1].xaxis_date()
        axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
        xlabels = axes[1].get_xticklabels()
        #pl.setp(xlabels, rotation=20)

        pl.show(block=False)
        fig.savefig(OUTFILE)

print('Saving statistics')
matstats_valid = {}
for run in runList:
    matstats_valid[run] = np.zeros((len(wmo_valid_list['P_l']),4))
    indw = 0 
    for iwmo,wmo in enumerate(wmo_list):
        if ~(matstats['P_l'][run][iwmo,0]==0):
            matstats_valid[run][indw,:] = matstats['P_l'][run][iwmo,:]
            indw += 1
ListMETRICS = ['IntRMS','IntBias','SurfRMS','SurfBias']
columnNames = ListMETRICS
rowNames = wmo_valid_list['P_l']
for run in runList:
    outfiletable = args.outdir + '/' + 'stats_floatP_l' + run + '.txt'
    writetable(outfiletable,matstats_valid[run],rowNames,columnNames)
    
