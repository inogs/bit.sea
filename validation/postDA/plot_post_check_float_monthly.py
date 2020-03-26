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
from commons import genUserDateList as DL
from commons import timerequestors as requestors
#from commons import Timelist
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


START_TIME = TL['P_l'].Timelist[0]
END___TIME = TL['P_l'].Timelist[-1]

MonthlyTL = DL.getTimeList(START_TIME,END___TIME,'months=1')
Ntot = TL['P_l'].nTimes
MonthList = TL['P_l'].getMonthlist()

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
    LISTmonthly = []
    arraymis = np.array(LISTchl[3+il])
    for mreq in MonthList:
        mind,_ = TL['P_l'].select(mreq)
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        LISTmonthly.append(np.nanmean(arraymis[maskm]))

    plt.plot(MonthlyTL,LISTmonthly,color=DICTcol[ll],label=ll)
    plt.plot(MonthlyTL[-1],LISTmonthly[-1],'o', \
             markeredgecolor='red',markerfacecolor=DICTcol[ll])

plt.legend(loc='upper left')
plt.grid()
plt.title(r'Misfit RMS chl $[mg/m^3]$')

plt.sca(axs[1])
LISTpercm = []
arrayobs = np.array(LISTchl[2])
arrayobsexc = np.array(LISTchl[5])
for mreq in MonthList:
    mind,_ = TL['P_l'].select(mreq)
    maskm = np.zeros(Ntot,dtype=np.bool)
    maskm[mind] = True
    percused = np.nansum(arrayobs[maskm])/  \
        (np.nansum(arrayobs[maskm])+np.nansum(arrayobsexc[maskm])) * 100.
    LISTpercm.append(percused)

plt.bar(MonthlyTL,LISTpercm,width=10,label='% Used obs ')
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
LISTprofN = []
LISTprofexcN = []
arrayprof = np.array(LISTchl[1])
arrayprofexc = np.array(LISTchl[6])
for mreq in MonthList:
    mind,_ = TL['N3n'].select(mreq)
    maskm = np.zeros(Ntot,dtype=np.bool)
    maskm[mind] = True
    nexcprofm = np.nansum(arrayprofexc[maskm])
    nprofm = np.nansum(arrayprof[maskm])
    LISTprofN.append(nprofm)
    LISTprofexcN.append(nexcprofm)

plt.bar(MonthlyTL,np.array(LISTprofN)+np.array(LISTprofexcN),\
    width=10,label='N profiles')
plt.bar(MonthlyTL[-1],LISTprofN[-1]+LISTprofexcN[-1],\
    width=10,color='red')
plt.bar(MonthlyTL,LISTprofexcN,width=10,label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'chl_monthly.png')


fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(NITlist):
    LISTmonthly = []
    arraymis = np.array(LISTnit[3+il])
    for mreq in MonthList:
        mind,_ = TL['N3n'].select(mreq)
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        LISTmonthly.append(np.nanmean(arraymis[maskm]))

    plt.plot(MonthlyTL,LISTmonthly,color=DICTcol[ll],label=ll)
    plt.plot(MonthlyTL[-1],LISTmonthly[-1],'o', \
             markeredgecolor='red',markerfacecolor=DICTcol[ll])

plt.legend(loc='upper left')
plt.grid()
plt.title(r'Misfit RMS nit $[mmol/m^3]$')

plt.sca(axs[1])
LISTpercm = []
for mreq in MonthList:
    mind,_ = TL['N3n'].select(mreq)
    maskm = np.zeros(Ntot,dtype=np.bool)
    maskm[mind] = True
    percused = np.nansum(np.array(LISTnit[2])[maskm])*1./  \
        (np.nansum(np.array(LISTnit[2])[maskm])+np.nansum(np.array(LISTnit[6])[maskm])) * 100.
    LISTpercm.append(percused)

plt.bar(MonthlyTL,LISTpercm,width=10,label='% Used obs ')
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
LISTprofN = []
LISTprofexcN = []
for mreq in MonthList:
    mind,_ = TL['N3n'].select(mreq)
    maskm = np.zeros(Ntot,dtype=np.bool)
    maskm[mind] = True
    nprofm = np.nansum(np.array(LISTnit[1])[maskm])
    nexcprofm = np.nansum(np.array(LISTnit[7])[maskm])
    LISTprofN.append(nprofm)
    LISTprofexcN.append(nexcprofm)

plt.bar(MonthlyTL,np.array(LISTprofN)+np.array(LISTprofexcN),\
    width=10,label='N profiles')
plt.bar(MonthlyTL[-1],LISTprofN[-1]+LISTprofexcN[-1],\
    width=10,color='red')
plt.bar(MonthlyTL,LISTprofexcN,width=10,label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'nit_monthly.png')



#plt.show(block=False)





