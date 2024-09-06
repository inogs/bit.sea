import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--indir', '-i',
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
from commons.utils import addsep

from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from instruments.matchup_manager import Matchup_Manager
from basins import OGS

from profiler_corr import *

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from timeseries.plot import *
import time
import pickle


INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)

RUN_LIST = ['HC_2017_simdd','HC_2017_assw']


font_s =  15 
label_s = 15


T_start = DATESTART
T_end   = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start,'%Y%m%d'))
T_end2num   = mpldates.date2num(datetime.strptime(T_end,'%Y%m%d'))


Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

print wmo_list

varRMSD_LIST = ['Chla','N3n','votemper','vosaline','density']

corrvarLIST = ['nitdens','nittemp']

DICTcolrun = {
    RUN_LIST[0]: 'orange',
    RUN_LIST[1]: 'green',
}

for j,wmo in enumerate(wmo_list):
# for j,wmo in enumerate([wmo_list[0]]):
    print(wmo)

    corrmod = {}
    covnmod = {}
    RMSD = {}
    for run in RUN_LIST:
        filename = INDIR + run + '/STATS/corrcov' + run + '_' + wmo + '.pkl'
        fid = open(filename,'r')
        LIST = pickle.load(fid)
        fid.close()
        corrmod[run] = LIST[1]
        covnmod[run] = LIST[2]
        RMSD[run] = LIST[3]

    filename = INDIR + run + '/STATS/corrcovFloat_' + wmo + '.pkl'
    fid = open(filename,'r')
    LISTF = pickle.load(fid)
    fid.close()

    timelist = LISTF[0]
    corrfloat = LISTF[1]
    covnfloat = LISTF[2]

    #Save txt
    timeymd = list()
    for tt in timelist:
        date8 = tt.strftime('%Y%m%d')
        timeymd.append(date8)
    ndates = len(timeymd)
    np.savetxt(OUTDIR + '/' + wmo + 'dates.txt',timeymd,fmt='%s')

    for run in RUN_LIST:

        RMSDarray = np.zeros((ndates,5))
        for ivar,var in enumerate(varRMSD_LIST):
            RMSDarray[:,ivar] = RMSD[run][var]
        np.savetxt(OUTDIR + '/' + wmo + 'RMSD_' + run + '.txt',RMSDarray)

        corrarray = np.zeros((ndates,2))
        for icorr,corrvar in enumerate(['nittemp','nitdens']):
            corrarray[:,icorr] = corrmod[run][corrvar]
        np.savetxt(OUTDIR + '/' + wmo + 'corr_' + run + '.txt',RMSDarray)

    corrarrayf = np.zeros((ndates,2))
    for icorr,corrvar in enumerate(['nittemp','nitdens']):
        corrarrayf[:,icorr] = corrfloat[corrvar]
    np.savetxt(OUTDIR + '/' + wmo + 'corr_float.txt',RMSDarray)


#     # Figures
#     plt.close('all')
#     fig,axes = plt.subplots(2,1,sharex=True)

#     plt.sca(axes[0])
#     plt.plot(timelist,corrfloat['nittemp'],label='Float')
#     for run in RUN_LIST:
#         plt.plot(timelist,corrmod[run]['nittemp'],label=run)
#     plt.legend(loc='best')
#     plt.grid()
#     plt.xlim(T_start2num,T_end2num)
#     plt.title('Corrcoef nit-temp ' + wmo)

#     plt.sca(axes[1])
#     plt.plot(timelist,corrfloat['nitdens'],label='Float')
#     for run in RUN_LIST:
#         plt.plot(timelist,corrmod[run]['nitdens'],label=run)
#     plt.legend(loc='best')
#     plt.grid()
#     plt.xlim(T_start2num,T_end2num)
#     plt.title('Corrcoef nit-dens ' + wmo)

#     fig.autofmt_xdate()

#     plt.savefig(OUTDIR + '/corrcoef' + wmo + '.png')


#     fig,axes = plt.subplots(2,1,sharex=True)

#     plt.sca(axes[0])
#     plt.plot(timelist,covnfloat['nittemp'],label='Float')
#     for run in RUN_LIST:
#         plt.plot(timelist,covnmod[run]['nittemp'],label=run)
#     plt.legend(loc='best')
#     plt.grid()
#     plt.xlim(T_start2num,T_end2num)
#     plt.title('Normalized covariance nit-temp' + wmo)

#     plt.sca(axes[1])
#     plt.plot(timelist,covnfloat['nitdens'],label='Float')
#     for run in RUN_LIST:
#         plt.plot(timelist,covnmod[run]['nitdens'],label=run)
#     plt.legend(loc='best')
#     plt.grid()
#     plt.xlim(T_start2num,T_end2num)
#     plt.title('Normalized covariance nit-dens ' + wmo)

#     fig.autofmt_xdate()

#     plt.savefig(OUTDIR + '/covnorm' + wmo + '.png')


#     fig,axes = plt.subplots(5,1,sharex=True,figsize=[8,8])

#     for ivar,var in enumerate(varRMSD_LIST):
#         plt.sca(axes[ivar])
#         for run in RUN_LIST:
#             plt.plot(timelist,RMSD[run][var],label=run,color=DICTcolrun[run])
#         plt.legend(loc='best')
#         plt.grid()
#         plt.xlim(T_start2num,T_end2num)
#         plt.title('RMSD ' + var + ' ' + wmo)

#     fig.autofmt_xdate()

#     plt.savefig(OUTDIR + '/RMSD' + wmo + '.png')


# #    plt.show(block=False)


