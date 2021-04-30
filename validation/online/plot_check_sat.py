import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    plot post check and mis stats of BGC-Argo DA
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--indir', '-i',
                        type=str,
                        default=None,
                        required=True,
                        help="Dir with output of post_check.py")

    parser.add_argument('--statsdir', '-s',
                        type=str,
                        default=None,
                        required=True,
                        help="Dir with output of sat_check")

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help="")

    parser.add_argument('--submaskdir', '-b',
                        type=str,
                        default=None,
                        required=True,
                        help="")

    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
import scipy.io.netcdf as NC
import commons.timerequestors as requestors
from commons.utils import addsep
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from basins import V2


INDIR = addsep(args.indir)
INCHECKDIR = addsep(args.statsdir)
OUTDIR = addsep(args.outdir)
SUBMASKDIR = addsep(args.submaskdir)

V2.P.basin_list = V2.P.basin_list[:-1] # removing Atlantic

TLmis = TimeList.fromfilenames(None,INDIR,'*.npy',
            prefix='',dateformat='%Y%m%d')
TI = TimeInterval(TLmis.Timelist[0].strftime('%Y%m%d'), \
                  TLmis.Timelist[-1].strftime('%Y%m%d'))

TLstats = TimeList.fromfilenames(TI, INCHECKDIR, \
          '*_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-*-v02_statsCHECK.nc', \
          prefix='',dateformat='%Y%m%d')


LISTmis = {}
LISTst = {}
Npsub = {}


for sub in V2.P:
    LISTmis[sub.name] = [[] for ll in range(4)]
    LISTst[sub.name] = [[] for ll in range(3)]
    filemasksub = SUBMASKDIR + 'masksub.' + sub.name + 'All.npy'
    masksub = np.load(filemasksub)
    Npsub[sub.name] = np.sum(masksub)


for filemis in TLmis.filelist:
    misfmean = np.load(filemis)
    for isub,sub in enumerate(V2.P):
        LISTmis[sub.name][0].append(misfmean[0,isub])
        LISTmis[sub.name][1].append(misfmean[1,isub])
        LISTmis[sub.name][2].append(misfmean[2,isub])
        LISTmis[sub.name][3].append(misfmean[3,isub])

for statsfile in TLstats.filelist:
    FS = NC.netcdf_file(statsfile,'r')
    vv = FS.variables['SUBstatistics_day'].data.copy()
    FS.close()
    for isub,sub in enumerate(V2.P):
        LISTst[sub.name][0].append(vv[isub,5])
        LISTst[sub.name][1].append(vv[isub,6]+vv[isub,7])
        LISTst[sub.name][2].append(Npsub[sub.name])

DICTcol = {
    'All': 'green',
    'Open': 'orange',
}

for sub in V2.P:
    #print sub.name
    plt.close('all')
    fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

    plt.sca(axs[0])
    for ii,masktype in enumerate(['All','Open']):
        plt.plot(TLmis.Timelist,LISTmis[sub.name][ii], \
                 color=DICTcol[masktype],label=masktype)
        plt.plot(TLmis.Timelist[-1],LISTmis[sub.name][ii][-1], \
                 'o',markeredgecolor='red',markerfacecolor=DICTcol[masktype])
    plt.legend(loc='upper left')
    plt.title(sub.name + r' Mean RMSD $[mgchl/m^3]$ - Last date '+ \
              TLmis.Timelist[-1].strftime('%Y%m%d'))
    plt.grid()

    plt.sca(axs[1])
    plt.bar(TLmis.Timelist,100.*np.array(LISTmis[sub.name][3]),width=3, \
            label='% Weekly coverage ' + \
                  '%1.f' %(100.*np.nanmean(LISTmis[sub.name][3])))
    plt.legend(loc='upper left')
    plt.grid()

    plt.sca(axs[2])
    plt.bar(TLstats.Timelist,LISTst[sub.name][2],width=4,label='Np subbasin')
    plt.bar(TLstats.Timelist,LISTst[sub.name][0],width=1,label='Obs')
    plt.bar(TLstats.Timelist,LISTst[sub.name][1],width=1,label='Excluded Obs') 
    plt.legend(loc='upper left')
    plt.grid()

    fig.autofmt_xdate()
    plt.savefig(OUTDIR + sub.name + 'misfmean.png')




