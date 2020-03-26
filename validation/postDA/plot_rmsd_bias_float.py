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
from layerinfo import DICTlayers


INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)

varLIST = ['P_l','N3n']

TL = {}
for var in varLIST:
    TL[var] = TimeList.fromfilenames(None, INDIR, 'poststats.*' + var + '*npy', \
            prefix='postcheck.',dateformat='%Y%m%d')

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

plt.close('all')


# RMSD figure
#chl
fig,axs = plt.subplots(Nlayers_chl,1,sharex=True,sharey=True,figsize=[10,9])
plt.suptitle(r'RMSD chl $[mg/m^3]$')

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
plt.suptitle(r'RMSD nit $[mmol/m^3]$')

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

plt.show(block=False)




