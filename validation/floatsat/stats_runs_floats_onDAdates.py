
import numpy as np
import matplotlib.pyplot as pl
import argparse
#import pickle
import matplotlib.dates as mdates

from basins import OGS
from instruments import lovbio_float as bio_float
from profileruns import runList,colorList
from profiler_floatsat import ALL_PROFILES,TL,T_INT
from profiler_floatsat import INPUTDIR as inputdir_run
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
from commons.Timelist import TimeList
from commons.utils import writetable
from commons.mask import Mask
from commons.dataextractor import DataExtractor

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

    parser.add_argument(   '--baseruns', '-R',
                            type = str,
                            required =True,
                            default = "/pico/scratch/userexternal/ateruzzi/",
                            help = ''' Output image dir'''
                            )
    
    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required = True,
                            help = '''maskfile (meshmask.nc)''')

    
    # parser.add_argument(   '--open_sea_file', '-o',
    #                         type = str,
    #                         required = True,
    #                         help = 'Input pickle file for open sea')

    return parser.parse_args()

args = argument()

TheMask=Mask(args.maskfile)
indz200 = TheMask.getDepthIndex(200) + 1

MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)
#############

DA_DIR = inputdir_run + '/../DA__FREQ_1/'
TL_weekly = TimeList.fromfilenames(T_INT,DA_DIR,'chl.*.nc',prefix='chl.')

matstats = {}
for run in runList:
    matstats[run] = np.zeros((len(wmo_list),4))

wmo_valid_list = []

for iwmo,wmo in enumerate(wmo_list):
# for iwmo,wmo in enumerate(['6901860']):
    print('Float ' + wmo)
    A = {}
    print('Reading input ----')
    wmo_track_list = bio_float.filter_by_wmo(ALL_PROFILES,wmo)
    nP = len(wmo_track_list)
    times = [p.time for p in wmo_track_list]
    Lons = [p.lon for p in wmo_track_list]
    Lats = [p.lat for p in wmo_track_list]
    nTime = len(times)

    for irun,run in enumerate(runList):
        print(run)
        INDIR = args.inputdir + '/RUN_' + run + '/tmp_nc/'
        INPUT_FILE = INDIR + wmo + ".nc"
        A[run] = ncreader(INPUT_FILE)
        _, ref = A[run].plotdata('P_l','Int_0-200')
        _, surf_ref = A[run].plotdata('P_l','SurfVal')
        model = np.zeros(nTime)
        surf_model = np.zeros(nTime)
        AVE_DIR = args.baseruns + args.inputdir + '/RUN_' + run + \
                    '/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/'
        TL_daily = TimeList.fromfilenames(T_INT,AVE_DIR, \
                                  'ave.*P_l.nc',prefix='ave.')
        for iit,tt in enumerate(times):
            tuesday_index = TL_weekly.find(tt)
            nearest_tuesday = TL_weekly.Timelist[tuesday_index]
            daily_index = TL_daily.find(nearest_tuesday)+1
            print(tt,daily_index)
            daily_file = TL_daily.filelist[daily_index]
            modval = DataExtractor(TheMask,daily_file,'P_l',dimvar=3).values
            lonfloat = wmo_track_list[iit].lon
            latfloat = wmo_track_list[iit].lat
            indi,indj = TheMask.convert_lon_lat_to_indices(lonfloat,latfloat)
            modprofile = modval[:indz200,indj,indi]
            surf_model[iit] = modprofile[0]
            model[iit] = np.nansum(modprofile*TheMask.dz[:indz200])/ \
                         TheMask.dz[:indz200].sum()

            
            
        if (~np.isnan(model).all() == True) & (~np.isnan(ref).all() == True):
            # pl.close('all')
            # fig,axes = pl.subplots(2,1,sharex=True,dpi=150)
            # fig.suptitle('P_l - ' + wmo,fontsize=12,color='b')
            # OUTFILE = args.outdir + "P_l_" + wmo + "DAdates.png"            
            matstats[run][iwmo,0] = (np.nanmean((model-ref)**2))**0.5
            matstats[run][iwmo,1] = np.nanmean(model-ref)
            matstats[run][iwmo,2] = (np.nanmean((surf_model-surf_ref)**2))**0.5
            matstats[run][iwmo,3] = np.nanmean(surf_model-surf_ref)
            if irun==0:
                wmo_valid_list.append(wmo)
        #         axes[0].plot(times,ref,'o-g',label='Float')
        #     axes[0].plot(times,model,'-',color=colorList[run],label=run)
        #     axes[0].set_title('Integral 0-200m',fontsize=10)

        #     if irun==0:
        #         axes[1].plot(times,surf_ref,'o-g',label='Float')
        #     axes[1].plot(times,surf_model,'-',color=colorList[run],label=run)
        #     axes[1].set_title('Surface',fontsize=10)
        #     axes[0].set_ylabel('Chlorophyll \n $[mg{\  } m^{-3}]$',fontsize=10)
        #     axes[1].set_ylabel('Chlorophyll \n $[mg{\  } m^{-3}]$',fontsize=10)
            
        # for iax in range(2):
        #     axes[iax].grid(True)
        #     legend0 = axes[iax].legend(loc='best', fontsize=8)
        #     axes[iax].tick_params(axis='both',labelsize=8)
        # axes[1].xaxis_date()
        # axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
        # xlabels = axes[1].get_xticklabels()
        #pl.setp(xlabels, rotation=20)

        #pl.show(block=False)
        #fig.savefig(OUTFILE)

# print('Saving statistics')
matstats_valid = {}
for run in runList:
    matstats_valid[run] = np.zeros((len(wmo_valid_list),4))
    indw = 0 
    for iwmo,wmo in enumerate(wmo_list):
        if ~(matstats[run][iwmo,0]==0):
            matstats_valid[run][indw,:] = matstats[run][iwmo,:]
            indw += 1
ListMETRICS = ['IntRMS','IntBias','SurfRMS','SurfBias']
columnNames = ListMETRICS
rowNames = wmo_valid_list
for run in runList:
    outfiletable = args.outdir + '/' + 'stats_floatP_l' + run + '_DAdates.txt'
    writetable(outfiletable,matstats_valid[run],rowNames,columnNames)
    
