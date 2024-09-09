import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help = "Directory where metrics are saved")

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--idwmo', '-w',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
import scipy.io.netcdf as NC

#from commons.layer import Layer
from commons.mask import Mask
from commons.utils import addsep
from metrics import *

from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from instruments.matchup_manager import Matchup_Manager

from profiler_compprofs import *

#from layer_integral import coastline
from basins import OGS

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from timeseries.plot import *
#import numpy.ma as ma
from mhelpers.pgmean import PLGaussianMean

from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader

TheMask = Mask(args.maskfile)

OUTDIR = addsep(args.outdir)

idwmo = args.idwmo

DIRSTATS = {}
for run in [RUN_REF,RUN_DA]:
	DIRSTATS[run] = addsep(args.indir) + '/' + run + '/STATS/'


font_s =  15 
label_s = 15


def readModelProfile(filename,var, wmo):
    ncIN = NC.netcdf_file(filename,'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile,:]
    return Profile

def get_level300(TheMask):

    i0 = TheMask.getDepthIndex(300)
    i1 = i0 + 1
      
    diff0 = 300 - TheMask.zlevels[i0]
    diff1 = 300 - TheMask.zlevels[i1]
    data = [(i0,diff0),(i1,diff1)]
    ix,datamin = min(data, key=lambda t: t[1])

    return ix


T_start = DATESTART
T_end   = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start,'%Y%m%d'))
T_end2num   = mpldates.date2num(datetime.strptime(T_end,'%Y%m%d'))
reg1 = [OGS.med]
reg_sn = ['med']

max_depth = get_level300(TheMask)
#max_depthp = get_level300(TheMask_Phys)
#layerStats = Layer(0,200)

VARLIST = ['Chla','N3n','O2o']
VARLIST = ['N3n']
VARLIST = ['P_l','N3n']
Adj = {
	'P_l':  True,
	'Chla': True,
	'O2o':  False,
	'N3n':  True,
}

nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11,1.0)

Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

MM = {}
for run in [RUN_REF,RUN_DA]:
    MM[run] = Matchup_Manager(ALL_PROFILES,TL,BASEDIR[run])

print wmo_list

DICTstyle = {
    RUN_REF: ['-b'],
    RUN_DA:  ['-r'],
    'Float': ['-g'],
}

#for j,wmo in enumerate(wmo_list):
for j,wmo in enumerate(wmo_list):
    if not(wmo==idwmo): continue
    print(wmo)

    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    nP = len(list_float_track)
    #NewPres_5m=np.linspace(0,300,61)
    #depths=NewPres_5m

    for ivar, var_mod in enumerate(VARLIST):
        print var_mod 
#        plotmat_model = {}
       # for run in [RUN_REF,RUN_DA]:
#            plotmat_model[run] = np.zeros([len(depths), len(list_float_track)])*np.nan
        var = LOVFLOATVARS[var_mod]
        adj=Adj[var_mod]

   
        for ip, p in enumerate(list_float_track):
            if not(var in p.available_params): continue
            Pres,Profile,_ = p.read(var,adj)
            #timelabel_list.append(p.time)
            date8 = p.time.strftime('%Y%m%d')
            print date8

            #continue
            plt.close('all')
            print '... figure'
            plt.figure(figsize=[4,5])

# PLOT FOR THE MODEL
            for run in [RUN_REF,RUN_DA]:
                print '      ' + run
                TM = MM[run].modeltime(p)
                FILENAME = BASEDIR[run] + TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
                M = readModelProfile(FILENAME,var_mod,p.ID())
                print M[:3]
                plt.plot(M[:max_depth+1], \
                     TheMask.zlevels[:max_depth+1], \
                     DICTstyle[run][0],label=run)

            plt.plot(Profile,Pres,DICTstyle['Float'][0],label='Float')

            plt.grid()
            plt.ylim(300,0)
            plt.legend(loc='best')

            if (var_mod == 'Chla') or (var_mod == 'P_l'):
                plt.xlim([0,.5])
            #    ax1.set_xlabel("$[mgchl/m^3]$")
            if (var_mod == 'N3n'):
                plt.xlim([0,6])
            #    ax1.set_xlabel("$[mmolN/m^3]$")

            plt.title(var_mod + ' ' + date8 + ' - wmo ' + wmo)

            #ax1.set_ylabel("depth $[m]$",color = 'k',fontsize=font_s)
            #plt.show(block=False)
            print (''.join([OUTDIR,'profs_',p.name(), \
                        '_',date8,'_',var_mod,'.png']))
            plt.savefig(''.join([OUTDIR,'profs_',p.name(), \
                        '_',date8,'_',var_mod,'.png']))

            plt.close('all')



