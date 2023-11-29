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

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from commons.utils import addsep
from commons.Timelist import TimeList
from commons.Timelist import TimeInterval
from commons import genUserDateList as DL
from layerinfo import DICTlayersQ

INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)

varLIST = ['P_l','N3n','O2o']

TL = {}
for var in varLIST:
    TL[var] = TimeList.fromfilenames(None, INDIR, 'postcheck.*' + var + '*npy', \
            prefix='postcheck.',dateformat='%Y%m%d')

START_TIME = TL['P_l'].Timelist[0]
END___TIME = TL['P_l'].Timelist[-1]

TI = TimeInterval(START_TIME.strftime('%Y%m%d'),END___TIME.strftime('%Y%m%d'))

MonthlyTL = DL.getTimeList(START_TIME,END___TIME,months=1)
Ntot = TL['P_l'].nTimes

mmTL = TimeList(MonthlyTL)
mmTL.inputFrequency='monthly' # added for tests with a single month
MonthList = mmTL.getOwnList()


LISTchl = [[] for ii in range(7)]
for dd in TL['P_l'].filelist:
    ll = np.load(dd, allow_pickle=True,encoding="latin1")
    for ii in range(7):
        LISTchl[ii].append(ll[ii])
    
LISTnit = [[] for ii in range(8)]
for dd in TL['N3n'].filelist:
    ll = np.load(dd, allow_pickle=True,encoding="latin1")
    for ii in range(8):
        LISTnit[ii].append(ll[ii])
    
LISToxy = [[] for ii in range(8)]
for dd in TL['O2o'].filelist:
    ll = np.load(dd, allow_pickle=True,encoding="latin1")
    for ii in range(8):
        LISToxy[ii].append(ll[ii])


CHLlist = ['0-50','50-150']
NITlist = ['0-50','50-150','300-400']
OXYlist = ['0-50','50-150','300-400']
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

##daily
# chl
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


# nit
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


# oxy
fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(OXYlist):
    plt.plot(LISToxy[0],LISToxy[3+il],label=ll,color=DICTcol[ll])
    plt.plot(LISToxy[0][-1],LISToxy[3+il][-1],'o', \
            markeredgecolor='red',markerfacecolor=DICTcol[ll])
plt.legend(loc='upper left')
plt.grid()
plt.title(r'Misfit RMS Oxy $[mmol/m^3]$ - Last date ' + \
        LISToxy[0][-1].strftime('%Y-%m-%d'))

plt.sca(axs[1])
percused = np.array(LISToxy[2])/ \
            (np.array(LISToxy[2])+np.array(LISToxy[6])) * 100.
percTOT = np.nansum(LISToxy[2])/ \
            (np.nansum(LISToxy[2]) + np.nansum(LISToxy[6])) *100.
plt.bar(LISToxy[0],percused,width=1,label='% Used obs ' + '%.1f' %percTOT)
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
plt.bar(LISToxy[0],np.array(LISToxy[1])+np.array(LISToxy[7]), \
    width=1,label='N profiles')
plt.bar(LISToxy[0],LISToxy[7],label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'oxy_daily.png')




## monthly
# chl
fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(CHLlist):
    LISTmonthly = []
    arraymis = np.array(LISTchl[3+il])
    for mreq in MonthList:
        mind,_ = TL['P_l'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntot,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arraymis[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arraymis[maskm]))
        else:
            LISTmonthly.append(np.nan)

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
    if len(mind)>0:
        maskm = np.zeros(Ntot,dtype=bool)
        maskm[mind] = True
        if np.all(np.isnan(arrayobs[maskm])):
            LISTpercm.append(np.nan)
        else:
            percused = np.nansum(arrayobs[maskm])/  \
                (np.nansum(arrayobs[maskm])+np.nansum(arrayobsexc[maskm])) * 100.
            LISTpercm.append(percused)
    else:
        LISTpercm.append(np.nan)

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
    if len(mind)>0:
        maskm = np.zeros(Ntot,dtype=bool)
        maskm[mind] = True
        nexcprofm = np.nansum(arrayprofexc[maskm])
        nprofm = np.nansum(arrayprof[maskm])
        LISTprofN.append(nprofm)
        LISTprofexcN.append(nexcprofm)
    else:
        LISTprofexcN.append(np.nan)
        LISTprofN.append(np.nan)

plt.bar(MonthlyTL,np.array(LISTprofN)+np.array(LISTprofexcN),\
    width=10,label='N profiles')
plt.bar(MonthlyTL[-1],LISTprofN[-1]+LISTprofexcN[-1],\
    width=10,color='red')
plt.bar(MonthlyTL,LISTprofexcN,width=10,label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'chl_monthly.png')
plt.close(fig)


#Nit
fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(NITlist):
    LISTmonthly = []
    arraymis = np.array(LISTnit[3+il])
    for mreq in MonthList:
        mind,_ = TL['N3n'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntot,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arraymis[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arraymis[maskm]))
        else :
            LISTmonthly.append(np.nan)

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
    if len(mind)>0:
        maskm = np.zeros(Ntot,dtype=bool)
        maskm[mind] = True
        if np.all(np.isnan(np.array(LISTnit[2])[maskm])):
            LISTpercm.append(np.nan)
        else:
            percused = np.nansum(np.array(LISTnit[2])[maskm])*1./  \
                (np.nansum(np.array(LISTnit[2])[maskm])+np.nansum(np.array(LISTnit[6])[maskm])) * 100.
            LISTpercm.append(percused)
    else :
        LISTpercm.append(np.nan)

plt.bar(MonthlyTL,LISTpercm,width=10,label='% Used obs ')
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
LISTprofN = []
LISTprofexcN = []

for mreq in MonthList:
    mind,_ = TL['N3n'].select(mreq)
    if len(mind)>0:
        maskm = np.zeros(Ntot,dtype=bool)
        maskm[mind] = True
        nprofm = np.nansum(np.array(LISTnit[1])[maskm])
        nexcprofm = np.nansum(np.array(LISTnit[7])[maskm])
        LISTprofN.append(nprofm)
        LISTprofexcN.append(nexcprofm)
    else :
        LISTprofexcN.append(np.nan)
        LISTprofN.append(np.nan)

plt.bar(MonthlyTL,np.array(LISTprofN)+np.array(LISTprofexcN),\
    width=10,label='N profiles')
plt.bar(MonthlyTL[-1],LISTprofN[-1]+LISTprofexcN[-1],\
    width=10,color='red')
plt.bar(MonthlyTL,LISTprofexcN,width=10,label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'nit_monthly.png')
plt.close(fig)

#oxy
fig,axs = plt.subplots(3,1,sharex=True,figsize=[10,6])

plt.sca(axs[0])
for il,ll in enumerate(OXYlist):
    LISTmonthly = []
    arraymis = np.array(LISToxy[3+il])
    for mreq in MonthList:
        mind,_ = TL['N3n'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntot,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arraymis[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arraymis[maskm]))
        else :
            LISTmonthly.append(np.nan)

    plt.plot(MonthlyTL,LISTmonthly,color=DICTcol[ll],label=ll)
    plt.plot(MonthlyTL[-1],LISTmonthly[-1],'o', \
             markeredgecolor='red',markerfacecolor=DICTcol[ll])

plt.legend(loc='upper left')
plt.grid()
plt.title(r'Misfit RMS oxy $[mmol/m^3]$')

plt.sca(axs[1])
LISTpercm = []
for mreq in MonthList:
    mind,_ = TL['N3n'].select(mreq)
    if len(mind)>0:
        maskm = np.zeros(Ntot,dtype=bool)
        maskm[mind] = True
        if np.all(np.isnan(np.array(LISToxy[2])[maskm])):
            LISTpercm.append(np.nan)
        else:
            percused = np.nansum(np.array(LISToxy[2])[maskm])*1./  \
                (np.nansum(np.array(LISToxy[2])[maskm])+np.nansum(np.array(LISToxy[6])[maskm])) * 100.
            LISTpercm.append(percused)
    else :
        LISTpercm.append(np.nan)

plt.bar(MonthlyTL,LISTpercm,width=10,label='% Used obs ')
plt.legend(loc='upper left')
plt.grid()

plt.sca(axs[2])
LISTprofN = []
LISTprofexcN = []

for mreq in MonthList:
    mind,_ = TL['N3n'].select(mreq)
    if len(mind)>0:
        maskm = np.zeros(Ntot,dtype=bool)
        maskm[mind] = True
        nprofm = np.nansum(np.array(LISToxy[1])[maskm])
        nexcprofm = np.nansum(np.array(LISToxy[7])[maskm])
        LISTprofN.append(nprofm)
        LISTprofexcN.append(nexcprofm)
    else :
        LISTprofexcN.append(np.nan)
        LISTprofN.append(np.nan)

plt.bar(MonthlyTL,np.array(LISTprofN)+np.array(LISTprofexcN),\
    width=10,label='N profiles')
plt.bar(MonthlyTL[-1],LISTprofN[-1]+LISTprofexcN[-1],\
    width=10,color='red')
plt.bar(MonthlyTL,LISTprofexcN,width=10,label='N profiles with exclusion')
plt.legend(loc='upper left')
plt.grid()

fig.autofmt_xdate()

plt.savefig(OUTDIR + 'oxy_monthly.png')
plt.close(fig)

## Quid layers
TLq = {}
for var in varLIST:
    TLq[var] = TimeList.fromfilenames(TI, INDIR, 'poststats.*' + var + '*npy', \
            prefix='poststats.',dateformat='%Y%m%d')

Ntotq = TLq['P_l'].nTimes
Nmonths = len(MonthlyTL)
# MonthList = TL['P_l'].getMonthlist()

CHLlist = DICTlayersQ['chl']
NITlist = DICTlayersQ['nit']
OXYlist = DICTlayersQ['oxy']
Nlayers_chl = len(CHLlist)
Nlayers_nit = len(NITlist)
Nlayers_oxy = len(OXYlist)

RMSDchl = [[] for ii in range(1+Nlayers_chl)]
BIASchl = [[] for ii in range(1+Nlayers_chl)]
for dd in TLq['P_l'].filelist:
    ll = np.load(dd, allow_pickle=True,encoding="latin1")
    RMSDchl[0].append(ll[0])
    BIASchl[0].append(ll[0])
    for ii in range(Nlayers_chl):
        RMSDchl[1+ii].append(ll[1+ii])
        BIASchl[1+ii].append(ll[1+ii+Nlayers_chl])
Nchl = len(RMSDchl[0])

RMSDnit = [[] for ii in range(1+Nlayers_nit)]
BIASnit = [[] for ii in range(1+Nlayers_nit)]
for dd in TLq['N3n'].filelist:
    ll = np.load(dd, allow_pickle=True, encoding="latin1")
    RMSDnit[0].append(ll[0])
    BIASnit[0].append(ll[0])
    for ii in range(Nlayers_nit):
        RMSDnit[1+ii].append(ll[1+ii])
        BIASnit[1+ii].append(ll[1+ii+Nlayers_nit])
Nnit = len(RMSDnit[0])

RMSDoxy = [[] for ii in range(1+Nlayers_oxy)]
BIASoxy = [[] for ii in range(1+Nlayers_oxy)]
for dd in TLq['O2o'].filelist:
    ll = np.load(dd, allow_pickle=True, encoding="latin1")
    RMSDoxy[0].append(ll[0])
    BIASoxy[0].append(ll[0])
    for ii in range(Nlayers_oxy):
        RMSDoxy[1+ii].append(ll[1+ii])
        BIASoxy[1+ii].append(ll[1+ii+Nlayers_oxy])
Noxy = len(RMSDoxy[0])

rangermsdCHL = {
    '0-10': [.04,.10,.06],
    '10-30': [.04,.10,.06],
    '30-60': [.05,.12,.07],
    '60-100': [.05,.07,.06],
    '100-150': [.04,.05,.05],
}

rangermsdNIT = {
    '0-10': [.21,.58,.46],
    '10-30': [.17,.59,.43],
    '30-60': [.32,.76,.47],
    '60-100': [.41,.65,.51],
    '100-150': [.34,.61,.43],
    '150-300': [.21,.57,.31],
    '300-600': [.11,.53,.28],
    '600-1000': [.29,1.17,.69],
}

rangermsdOXY = {
    '0-10': [2.1,6.8,3.6],
    '10-30': [2.9,9.1,4.7],
    '30-60': [5.1,10.1,6.6],
    '60-100': [4.1,8.3,6.2],
    '100-150': [4,7.4,5.7],
    '150-300': [3.1,7.5,4.6],
    '300-600': [0.8,6.9,4.0],
    '600-1000': [2.7,8.8,5.9],
}

rangebiasCHL = {
    '0-10': [-.06,.0,-.03],
    '10-30': [-.06,0.0,-.03],
    '30-60': [-.08,-.01,-.04],
    '60-100': [-.06,.01,-.01],
    '100-150': [.01,.03,0.025],
}

rangebiasNIT = {
    '0-10': [-.16,.50,.15],
    '10-30': [-.20,.43,.12],
    '30-60': [-.36,.29,-.03],
    '60-100': [-.36,.40,-.06],
    '100-150': [-.18,.43,.01],
    '150-300': [-.13,.11,.01],
    '300-600': [-.36,-.07,-.17],
    '600-1000': [-1.15,-.29,-.74],
}

rangebiasOXY = {
    '0-10': [-2.6,.7,-.91],
    '10-30': [-5.1,-1,-2.57],
    '30-60': [-6,-0.1,-2.64],
    '60-100': [-1.8,5.2,-1.17],
    '100-150': [-3,2.7,-.53],
    '150-300': [-6.7,3.1,-.38],
    '300-600': [-.1,5.9,3.13],
    '600-1000': [-2.6,8,4.21],
}


# Daily
plt.close('all')


# RMSD figure
#chl
fig,axs = plt.subplots(Nlayers_chl,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'Misfit RMS chl $[mg/m^3]$')

for il,ll in enumerate(CHLlist):
    plt.sca(axs[il])
    #for im in range(2):
    #    plt.plot(RMSDchl[0],[rangermsdCHL[ll][im] for ii in range(Nchl)], \
    #             ':',color='grey')
    plt.plot(RMSDchl[0],[rangermsdCHL[ll][2] for ii in range(Nchl)], \
             '-',color='grey',label='QuID mean %.3f' %rangermsdCHL[ll][2])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDchl[1+il]))
    plt.plot(RMSDchl[0],RMSDchl[1+il],label=txtlabel)
    plt.plot(RMSDchl[0][-1],RMSDchl[1+il][-1],'o',color='red')
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'chl_rmsdlayers.png')

#nit
fig,axs = plt.subplots(Nlayers_nit,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'Misfit RMS nit $[mmol/m^3]$')

for il,ll in enumerate(NITlist):
    plt.sca(axs[il])
    #for im in range(2):
    #    plt.plot(RMSDnit[0],[rangermsdNIT[ll][im] for ii in range(Nnit)], \
    #             ':',color='grey')
    plt.plot(RMSDnit[0],[rangermsdNIT[ll][2] for ii in range(Nnit)], \
             '-',color='grey',label='QuID mean %.3f' %rangermsdNIT[ll][2])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDnit[1+il]))
    plt.plot(RMSDnit[0],RMSDnit[1+il],label=txtlabel)
    plt.plot(RMSDnit[0][-1],RMSDnit[1+il][-1],'o',color='red')
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'nit_rmsdlayers.png')

#oxy
fig,axs = plt.subplots(Nlayers_oxy,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'Misfit RMS oxy $[mmol/m^3]$')

for il,ll in enumerate(OXYlist):
    plt.sca(axs[il])
    #for im in range(2):
    #    plt.plot(RMSDoxy[0],[rangermsdoxy[ll][im] for ii in range(Noxy)], \
    #             ':',color='grey')
    plt.plot(RMSDoxy[0],[rangermsdOXY[ll][2] for ii in range(Noxy)], \
             '-',color='grey',label='QuID mean %.3f' %rangermsdOXY[ll][2])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDoxy[1+il]))
    plt.plot(RMSDoxy[0],RMSDoxy[1+il],label=txtlabel)
    plt.plot(RMSDoxy[0][-1],RMSDoxy[1+il][-1],'o',color='red')
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'oxy_rmsdlayers.png')


# BIAS figure
#chl
fig,axs = plt.subplots(Nlayers_chl,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'BIAS chl $[mg/m^3]$')

for il,ll in enumerate(CHLlist):
    plt.sca(axs[il])
    #for im in range(2):
    #    plt.plot(BIASchl[0],[rangebiasCHL[ll][im] for ii in range(Nchl)], \
    #             ':',color='grey')
    plt.plot(BIASchl[0],[rangebiasCHL[ll][2] for ii in range(Nchl)], \
             '-',color='grey',label='QuID mean %.3f' %rangebiasCHL[ll][2])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(BIASchl[1+il]))
    plt.plot(BIASchl[0],BIASchl[1+il],label=txtlabel)
    plt.plot(BIASchl[0][-1],BIASchl[1+il][-1],'o',color='red')
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'chl_biaslayers.png')


#nit
fig,axs = plt.subplots(Nlayers_nit,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'BIAS nit $[mmol/m^3]$')

for il,ll in enumerate(NITlist):
    plt.sca(axs[il])
    #for im in range(2):
    #    plt.plot(BIASnit[0],[rangebiasNIT[ll][im] for ii in range(Nnit)], \
    #             ':',color='grey')
    plt.plot(BIASnit[0],[rangebiasNIT[ll][2] for ii in range(Nnit)], \
             '-',color='grey',label='QuID mean %.3f' %rangebiasNIT[ll][2])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(BIASnit[1+il]))
    plt.plot(BIASnit[0],BIASnit[1+il],label=txtlabel)
    plt.plot(BIASnit[0][-1],BIASnit[1+il][-1],'o',color='red')
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'nit_biaslayers.png')

#oxy
fig,axs = plt.subplots(Nlayers_oxy,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'BIAS oxy $[mmol/m^3]$')

for il,ll in enumerate(OXYlist):
    plt.sca(axs[il])
    #for im in range(2):
    #    plt.plot(BIASoxy[0],[rangebiasoxy[ll][im] for ii in range(Noxy)], \
    #             ':',color='grey')
    plt.plot(BIASoxy[0],[rangebiasOXY[ll][2] for ii in range(Noxy)], \
             '-',color='grey',label='QuID mean %.3f' %rangebiasOXY[ll][2])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(BIASoxy[1+il]))
    plt.plot(BIASoxy[0],BIASoxy[1+il],label=txtlabel)
    plt.plot(BIASoxy[0][-1],BIASoxy[1+il][-1],'o',color='red')
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'oxy_biaslayers.png')

# Monthly

# RMSD figure
#chl
fig,axs = plt.subplots(Nlayers_chl,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'Misfit RMS chl $[mg/m^3]$')

for il,ll in enumerate(CHLlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDchl[1+il]))
    LISTmonthly = []
    arrayrms = np.array(RMSDchl[1+il])
    for mreq in MonthList:
        mind,_ = TLq['P_l'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntotq,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arrayrms[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arrayrms[maskm]))
        else:
            LISTmonthly.append(np.nan)

    plt.bar(MonthlyTL,LISTmonthly,width=10,label=txtlabel)
    plt.bar(MonthlyTL[-1],LISTmonthly[-1],width=10,color='red')
    #for im in range(2):
    #    plt.plot(MonthlyTL,[rangermsdCHL[ll][im] for ii in range(Nmonths)], \
    #             ':',color='grey')
    plt.plot(MonthlyTL,[rangermsdCHL[ll][2] for ii in range(Nmonths)], \
             '-',color='grey',label='QuID mean %.3f' %rangermsdCHL[ll][2])
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'chl_rmsdlayers_monthly.png')



#nit
fig,axs = plt.subplots(Nlayers_nit,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'Misfit RMS nit $[mmol/m^3]$')

for il,ll in enumerate(NITlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDnit[1+il]))
    LISTmonthly = []
    arrayrms = np.array(RMSDnit[1+il])
    for mreq in MonthList:
        mind,_ = TLq['N3n'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntotq,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arrayrms[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arrayrms[maskm]))
        else:
            LISTmonthly.append(np.nan)

    plt.bar(MonthlyTL,LISTmonthly,width=10,label=txtlabel)
    plt.bar(MonthlyTL[-1],LISTmonthly[-1],width=10,color='red')
    #for im in range(2):
    #    plt.plot(MonthlyTL,[rangermsdNIT[ll][im] for ii in range(Nmonths)], \
    #             ':',color='grey')
    plt.plot(MonthlyTL,[rangermsdNIT[ll][2] for ii in range(Nmonths)], \
             '-',color='grey',label='QuID mean %.3f' %rangermsdNIT[ll][2])
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'nit_rmsdlayers_monthly.png')

#oxy
fig,axs = plt.subplots(Nlayers_oxy,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'Misfit RMS oxy $[mmol/m^3]$')

for il,ll in enumerate(OXYlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDoxy[1+il]))
    LISTmonthly = []
    arrayrms = np.array(RMSDoxy[1+il])
    for mreq in MonthList:
        mind,_ = TLq['N3n'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntotq,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arrayrms[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arrayrms[maskm]))
        else:
            LISTmonthly.append(np.nan)

    plt.bar(MonthlyTL,LISTmonthly,width=10,label=txtlabel)
    plt.bar(MonthlyTL[-1],LISTmonthly[-1],width=10,color='red')
    #for im in range(2):
    #    plt.plot(MonthlyTL,[rangermsdoxy[ll][im] for ii in range(Nmonths)], \
    #             ':',color='grey')
    plt.plot(MonthlyTL,[rangermsdOXY[ll][2] for ii in range(Nmonths)], \
             '-',color='grey',label='QuID mean %.3f' %rangermsdOXY[ll][2])
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'oxy_rmsdlayers_monthly.png')


# BIAS figure
#chl
fig,axs = plt.subplots(Nlayers_chl,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'BIAS chl $[mg/m^3]$')

for il,ll in enumerate(CHLlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(BIASchl[1+il]))
    LISTmonthly = []
    arrayrms = np.array(BIASchl[1+il])
    for mreq in MonthList:
        mind,_ = TLq['P_l'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntotq,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arrayrms[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arrayrms[maskm]))
        else:
            LISTmonthly.append(np.nan)


    plt.bar(MonthlyTL,LISTmonthly,width=10,label=txtlabel)
    plt.bar(MonthlyTL[-1],LISTmonthly[-1],width=10,color='red')
    #for im in range(2):
    #    plt.plot(MonthlyTL,[rangebiasCHL[ll][im] for ii in range(Nmonths)], \
    #             ':',color='grey')
    plt.plot(MonthlyTL,[rangebiasCHL[ll][2] for ii in range(Nmonths)], \
             '-',color='grey',label='QuID mean %.3f' %rangebiasCHL[ll][2])
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'chl_biaslayers_monthly.png')



#nit
fig,axs = plt.subplots(Nlayers_nit,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'BIAS nit $[mmol/m^3]$')

for il,ll in enumerate(NITlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(BIASnit[1+il]))
    LISTmonthly = []
    arrayrms = np.array(BIASnit[1+il])
    for mreq in MonthList:
        mind,_ = TLq['N3n'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntotq,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arrayrms[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arrayrms[maskm]))
        else:
            LISTmonthly.append(np.nan)

    plt.bar(MonthlyTL,LISTmonthly,width=10,label=txtlabel)
    plt.bar(MonthlyTL[-1],LISTmonthly[-1],width=10,color='red')
    #for im in range(2):
    #    plt.plot(MonthlyTL,[rangebiasNIT[ll][im] for ii in range(Nmonths)], \
    #             ':',color='grey')
    plt.plot(MonthlyTL,[rangebiasNIT[ll][2] for ii in range(Nmonths)], \
             '-',color='grey',label='QuID mean %.3f' %rangebiasNIT[ll][2])
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'nit_biaslayers_monthly.png')


#oxy
fig,axs = plt.subplots(Nlayers_oxy,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'BIAS oxy $[mmol/m^3]$')

for il,ll in enumerate(OXYlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(BIASoxy[1+il]))
    LISTmonthly = []
    arrayrms = np.array(BIASoxy[1+il])
    for mreq in MonthList:
        mind,_ = TLq['N3n'].select(mreq)
        if len(mind)>0:
            maskm = np.zeros(Ntotq,dtype=bool)
            maskm[mind] = True
            if np.all(np.isnan(arrayrms[maskm])):
                LISTmonthly.append(np.nan)
            else:
                LISTmonthly.append(np.nanmean(arrayrms[maskm]))
        else:
            LISTmonthly.append(np.nan)

    plt.bar(MonthlyTL,LISTmonthly,width=10,label=txtlabel)
    plt.bar(MonthlyTL[-1],LISTmonthly[-1],width=10,color='red')
    #for im in range(2):
    #    plt.plot(MonthlyTL,[rangebiasoxy[ll][im] for ii in range(Nmonths)], \
    #             ':',color='grey')
    plt.plot(MonthlyTL,[rangebiasOXY[ll][2] for ii in range(Nmonths)], \
             '-',color='grey',label='QuID mean %.3f' %rangebiasOXY[ll][2])
    plt.legend(loc='upper left')
    plt.grid()

fig.autofmt_xdate()
plt.savefig(OUTDIR + 'oxy_biaslayers_monthly.png')



