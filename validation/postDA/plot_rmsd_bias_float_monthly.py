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
from commons import genUserDateList as DL
from layerinfo import DICTlayers


INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)

varLIST = ['P_l','N3n']

TL = {}
for var in varLIST:
    TL[var] = TimeList.fromfilenames(None, INDIR, 'poststats.*' + var + '*npy', \
            prefix='postcheck.',dateformat='%Y%m%d')

START_TIME = TL['P_l'].Timelist[0]
END___TIME = TL['P_l'].Timelist[-1]

MonthlyTL = DL.getTimeList(START_TIME,END___TIME,'months=1')
Ntot = TL['P_l'].nTimes
Nmonths = len(MonthlyTL)
MonthList = TL['P_l'].getMonthlist()



CHLlist = DICTlayers['chl']
NITlist = DICTlayers['nit']
Nlayers_chl = len(CHLlist)
Nlayers_nit = len(NITlist)

RMSDchl = [[] for ii in range(1+Nlayers_chl)]
BIASchl = [[] for ii in range(1+Nlayers_chl)]
for dd in TL['P_l'].filelist:
    ll = np.load(dd)
    RMSDchl[0].append(ll[0])
    BIASchl[0].append(ll[0])
    for ii in range(Nlayers_chl):
        RMSDchl[1+ii].append(ll[1+ii])
        BIASchl[1+ii].append(ll[1+ii+Nlayers_chl])
Nchl = len(RMSDchl[0])
    
RMSDnit = [[] for ii in range(1+Nlayers_nit)]
BIASnit = [[] for ii in range(1+Nlayers_nit)]
for dd in TL['N3n'].filelist:
    ll = np.load(dd)
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
    '60-100': [-.065,.018,-.024],
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

plt.close('all')


# RMSD figure
#chl
fig,axs = plt.subplots(Nlayers_chl,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'RMSD chl $[mg/m^3]$')

for il,ll in enumerate(CHLlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDchl[1+il]))
    LISTmonthly = []
    arrayrms = np.array(RMSDchl[1+il])
    for mreq in MonthList:
        mind,_ = TL['P_l'].select(mreq)
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        LISTmonthly.append(np.nanmean(arrayrms[maskm]))

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
plt.suptitle(r'RMSD nit $[mmol/m^3]$')

for il,ll in enumerate(NITlist):
    plt.sca(axs[il])
    txtlabel = ll +  ' - Mean %.3f' %(np.nanmean(RMSDnit[1+il]))
    LISTmonthly = []
    arrayrms = np.array(RMSDnit[1+il])
    for mreq in MonthList:
        mind,_ = TL['N3n'].select(mreq)
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        LISTmonthly.append(np.nanmean(arrayrms[maskm]))

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
        mind,_ = TL['P_l'].select(mreq)
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        LISTmonthly.append(np.nanmean(arrayrms[maskm]))

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
        mind,_ = TL['N3n'].select(mreq)
        maskm = np.zeros(Ntot,dtype=np.bool)
        maskm[mind] = True
        LISTmonthly.append(np.nanmean(arrayrms[maskm]))

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

plt.show(block=False)




