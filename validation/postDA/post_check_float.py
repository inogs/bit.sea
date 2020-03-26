import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    post of check and mis for BGC-Argo DA
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--incheckdir', '-i',
                        type=str,
                        default=None,
                        required=True,
                        help="Dir with output of check")

    parser.add_argument('--indadir', '-d',
                        type=str,
                        default=None,
                        required=True,
                        help="DA dir DA__FREQ_1")

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help="")

    return parser.parse_args()


args = argument()

import numpy as np
import commons.timerequestors as requestors
import os

from commons.utils import addsep
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList




INCHECKDIR = addsep(args.incheckdir)
INDADIR = addsep(args.indadir)
OUTDIR = addsep(args.outdir)



varLIST = ['chl','nit']
DICTvarname = {
    'chl': 'P_l',
    'nit': 'N3n',
}

TL = {}
for var in varLIST:
    TL[var] = TimeList.fromfilenames(None, INCHECKDIR, \
              '*' + DICTvarname[var] + '_check.txt', \
              prefix='',dateformat='%Y%m%d')

TLmis = TimeList.fromfilenames(None, INDADIR,  \
             '????????.arg_mis.dat', \
             prefix='',dateformat='%Y%m%d')



DICTlayers = {
    'chl': ['0-50','50-150'],
    'nit': ['0-50','50-150','300-400'],
}

DICTlayer = {
    '0-50': [0,50],
    '50-150': [50,150],
    '300-400': [300,400],
}
dep5m ={
    '0-50': range(0,51,5),
    '50-150': range(50,151,5),
    '300-400': range(300,401,5),
}

nexc = {}
Nlayers = {}
for var in varLIST:
    Nlayers[var] = len(DICTlayers[var])
    nexc[var] = 0

DICTflagvar = {
    'chl': 0,
    'nit': 1,
}

Ndates = TLmis.nTimes
for misfile,datemis in zip(TLmis.filelist,TLmis.Timelist):
    date8 = datemis.strftime('%Y%m%d')
    print misfile
    req = requestors.Daily_req(datemis.year,datemis.month,datemis.day)
    misALL = np.loadtxt(misfile,skiprows=1)
    LIST = {}
    for var in varLIST:
        LIST[var] = [i for i in range(5+Nlayers[var])]
        LIST[var][0] = datemis
        misvar = misALL[misALL[:,1]==DICTflagvar[var]]
    
        wmovar = set(misvar[:,-1])
        LIST[var][1]= len(wmovar)
        if len(misvar)>0:
            LIST[var][2] = len(misvar)
        else:
            LIST[var][2] = np.nan

        for il,ll in enumerate(DICTlayers[var]):
        # print '   ...' + ll
            LISTmeanmis = []
            for wmo in wmovar:
                # print wmo
                miswmo = misvar[misvar[:,-1]==wmo]
                maskll = (miswmo[:,4]>=DICTlayer[ll][0]) & \
                         (miswmo[:,4]<DICTlayer[ll][1])
                if any(maskll):
                    depmasked = miswmo[:,4][maskll]
                    mismasked = miswmo[:,6][maskll]
                    misint = np.interp(dep5m[ll],depmasked,mismasked)
                    rmsd = (np.nanmean(misint**2))**.5
                    LISTmeanmis.append(rmsd)
            if len(LISTmeanmis)>0:
                meanmis = np.nanmean(LISTmeanmis)
            else:
                meanmis = np.nan
            LIST[var][3+il] = meanmis


        icheck = TL[var].select_one(req)
        txtvar = TL[var].filelist[icheck]
        ss = os.path.getsize(txtvar)
        if ss==0:
            nobsexc = 0
            nprofexc = 0
        else:
            checkvar = np.loadtxt(txtvar,ndmin=2)
            nobsexc = np.sum(checkvar[:,5])
            nprofexc = checkvar.shape[0]
            nexc[var]+= 1 
            print ' ...some exclusions on ' + var + ' ' + np.str(nexc[var])
        LIST[var][3+Nlayers[var]] = nobsexc
        LIST[var][4+Nlayers[var]] = nprofexc

        nomefile = 'postcheck.' + date8 + '.' + DICTvarname[var] + '.npy'
        np.save(OUTDIR + nomefile, LIST[var])


print ' TOT exc CHL ' + np.str(nexc['chl']) + ' on ' + np.str(Ndates)
print ' TOT exc N3n ' + np.str(nexc['nit']) + ' on ' + np.str(Ndates)


# for var in varLIST:
    # nomefile = 'postcheck_' + var + '_'  + START_TIME + '_' + END___TIME
    # print nomefile
    # np.save(OUTDIR + nomefile,LIST[var])



