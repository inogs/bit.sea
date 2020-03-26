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

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help="")

    return parser.parse_args()

args = argument()

import numpy as np
import pylab as plt

from commons.utils import addsep
from commons.Timelist import TimeList


INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)

varLIST = ['P_l','N3n']

TL = {}
for var in varLIST:
    TL[var] = TimeList.fromfilenames(None, INDIR, 'postcheck.*' + var + '*npy', \
            prefix='postcheck.',dateformat='%Y%m%d')

LISTchl = [[] for ii in range(7)]
for dd in TL['P_l'].filelist:
    ll = np.load(dd)
    for ii in range(7):
        LISTchl[ii].append(ll[ii])
    
LISTnit = [[] for ii in range(8)]
for dd in TL['N3n'].filelist:
    ll = np.load(dd)
    for ii in range(8):
        LISTnit[ii].append(ll[ii])
    

CHLlist = ['0-50','50-150']
NITlist = ['0-50','50-150','300-400']
DICTlayer = {
    '0-50': [0,50],
    '50-150': [50,150],
    '300-400': [300,400],
}

DICTcol = {
    '0-50': 'green',
    '50-150': 'orange',
    '300-400': 'teal',
}


plt.close('all')

fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(CHLlist):
    plt.plot(LISTchl[0],LISTchl[3+il],label=ll,color=DICTcol[ll])
    plt.plot(LISTchl[0][-1],LISTchl[3+il][-1],'o', \
            markeredgecolor='red',markerfacecolor=DICTcol[ll])
plt.legend(loc='upper left')
plt.grid()
plt.title(r'Misfit RMS chl $[mg/m^3]$ - Last date '  + \
        LISTchl[0][-1].strftime('%Y-%m-%d'))

plt.sca(axs[1])
percused = np.array(LISTchl[2])/ \
           (np.array(LISTchl[2])+np.array(LISTchl[5])) * 100.
percTOT = np.nansum(LISTchl[2])/ \
           (np.nansum(LISTchl[2]) + np.nansum(LISTchl[5])) *100.
plt.bar(LISTchl[0],percused,width=1,label='% Used obs ' + '%.1f' %percTOT)
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
plt.bar(LISTchl[0],np.array(LISTchl[1])+np.array(LISTchl[6]),\
    width=1,label='N profiles')
plt.bar(LISTchl[0],LISTchl[6],label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'chl_daily.png')



fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(NITlist):
    plt.plot(LISTnit[0],LISTnit[3+il],label=ll,color=DICTcol[ll])
    plt.plot(LISTnit[0][-1],LISTnit[3+il][-1],'o', \
            markeredgecolor='red',markerfacecolor=DICTcol[ll])
plt.legend(loc='upper left')
plt.grid()
plt.title(r'Misfit RMS Nit $[mmol/m^3]$ - Last date ' + \
        LISTnit[0][-1].strftime('%Y-%m-%d'))

plt.sca(axs[1])
percused = np.array(LISTnit[2])/ \
            (np.array(LISTnit[2])+np.array(LISTnit[6])) * 100.
percTOT = np.nansum(LISTnit[2])/ \
            (np.nansum(LISTnit[2]) + np.nansum(LISTnit[6])) *100.
plt.bar(LISTnit[0],percused,width=1,label='% Used obs ' + '%.1f' %percTOT)
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
plt.bar(LISTnit[0],np.array(LISTnit[1])+np.array(LISTnit[7]), \
    width=1,label='N profiles')
plt.bar(LISTnit[0],LISTnit[7],label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'nit_daily.png')


plt.show(block=False)





