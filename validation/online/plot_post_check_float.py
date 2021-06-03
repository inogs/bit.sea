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

varLIST = ['P_l','N3n']

TL = {}
for var in varLIST:
    TL[var] = TimeList.fromfilenames(None, INDIR, 'postcheck.*' + var + '*npy', \
            prefix='postcheck.',dateformat='%Y%m%d')

START_TIME = TL['P_l'].Timelist[0]
END___TIME = TL['P_l'].Timelist[-1]

TI = TimeInterval(START_TIME.strftime('%Y%m%d'),END___TIME.strftime('%Y%m%d'))

MonthlyTL = DL.getTimeList(START_TIME,END___TIME,'months=1')
Ntot = TL['P_l'].nTimes

mmTL = TimeList(MonthlyTL)
MonthList = mmTL.getOwnList()


LISTchl = [[] for ii in range(7)]
for dd in TL['P_l'].filelist:
    ll = np.load(dd, allow_pickle=True)
    for ii in range(7):
        LISTchl[ii].append(ll[ii])
    
LISTnit = [[] for ii in range(8)]
for dd in TL['N3n'].filelist:
    ll = np.load(dd, allow_pickle=True)
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
            maskm = np.zeros(Ntot,dtype=np.bool)
            maskm[mind] = True
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
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        percused = np.nansum(arrayobs[maskm])/  \
            (np.nansum(arrayobs[maskm])+np.nansum(arrayobsexc[maskm])) * 100.
        LISTpercm.append(percused)
    else:
        LISTpercm.append(np.nanmean)

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
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        nexcprofm = np.nansum(arrayprofexc[maskm])
        nprofm = np.nansum(arrayprof[maskm])
        LISTprofN.append(nprofm)
        LISTprofexcN.append(nexcprofm)
    else:
        LISTprofexcN.append(np.nan)

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
            maskm = np.zeros(Ntot,dtype=np.bool)
            maskm[mind] = True
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
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
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
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        nprofm = np.nansum(np.array(LISTnit[1])[maskm])
        nexcprofm = np.nansum(np.array(LISTnit[7])[maskm])
        LISTprofN.append(nprofm)
        LISTprofexcN.append(nexcprofm)
    else :
        LISTprofexcN.append(np.nan)

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
Nlayers_chl = len(CHLlist)
Nlayers_nit = len(NITlist)

RMSDchl = [[] for ii in range(1+Nlayers_chl)]
BIASchl = [[] for ii in range(1+Nlayers_chl)]
for dd in TLq['P_l'].filelist:
    ll = np.load(dd, allow_pickle=True)
    RMSDchl[0].append(ll[0])
    BIASchl[0].append(ll[0])
    for ii in range(Nlayers_chl):
        RMSDchl[1+ii].append(ll[1+ii])
        BIASchl[1+ii].append(ll[1+ii+Nlayers_chl])
Nchl = len(RMSDchl[0])

RMSDnit = [[] for ii in range(1+Nlayers_nit)]
BIASnit = [[] for ii in range(1+Nlayers_nit)]
for dd in TLq['N3n'].filelist:
    ll = np.load(dd, allow_pickle=True)
    RMSDnit[0].append(ll[0])
    BIASnit[0].append(ll[0])
    for ii in range(Nlayers_nit):
        RMSDnit[1+ii].append(ll[1+ii])
        BIASnit[1+ii].append(ll[1+ii+Nlayers_nit])
Nnit = len(RMSDnit[0])

rangermsdCHL = {
    '0-10': [.039,.075,.055],
    '10-30': [.039,.076,.055],
    '30-60': [.039,.085,.063],
    '60-100': [.037,.086,.064],
    '100-150': [.03,.055,.041],
}

rangermsdNIT = {
    '0-10': [.15,.4,.29],
    '10-30': [.12,.37,.27],
    '30-60': [.2,.42,.3],
    '60-100': [.25,.59,.44],
    '100-150': [.41,.78,.56],
    '150-300': [.29,.44,.33],
    '300-600': [.53,1.08,.85],
    '600-1000': [.26,1.35,.9],
}

rangebiasCHL = {
    '0-10': [-.042,.015,-.016],
    '10-30': [-.042,0.019,-.0165],
    '30-60': [-.041,.014,-.0185],
    '60-100': [-.065,.018,.024],
    '100-150': [-.019,.011,-0.007],
}

rangebiasNIT = {
    '0-10': [-.19,.4,.14],
    '10-30': [-.21,.37,.12],
    '30-60': [-.34,.21,-.02],
    '60-100': [-.48,.32,-.016],
    '100-150': [-.25,.46,.17],
    '150-300': [-.18,.2,.032],
    '300-600': [-.91,.92,-.28],
    '600-1000': [-1.35,.21,-.78],
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
            maskm = np.zeros(Ntotq,dtype=np.bool)
            maskm[mind] = True
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
            maskm = np.zeros(Ntotq,dtype=np.bool)
            maskm[mind] = True
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
            maskm = np.zeros(Ntotq,dtype=np.bool)
            maskm[mind] = True
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
            maskm = np.zeros(Ntotq,dtype=np.bool)
            maskm[mind] = True
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




