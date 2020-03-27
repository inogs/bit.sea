import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    post of check and mis for BGC-Argo DA
    ''', formatter_class=argparse.RawTextHelpFormatter)

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
from layerinfo import *




INDADIR = addsep(args.indadir)
OUTDIR = addsep(args.outdir)



varLIST = ['chl','nit']
DICTvarname = {
    'chl': 'P_l',
    'nit': 'N3n',
}


TLmis = TimeList.fromfilenames(None, INDADIR,  \
             '????????.arg_mis.dat', \
             prefix='',dateformat='%Y%m%d')


#LayerLIST = [
#    [0,10],
#    [10,30],
#    [30,60],
#    [60,100],
#    [100,150],
#    [150,300],
#    [300,600]
#]

#DICTlayer = {}
#dep5m = {}
#DICTlayers = {
#    'chl' : [],
#    'nit' : [],
#}
#for ll in LayerLIST:
#    layername = '%s' %(ll[0]) + '-' + '%s' %(ll[1])
#    DICTlayer[layername] = ll
#    dep5m[layername] = range(ll[0],ll[1]+1,5)
#    DICTlayers['nit'].append(layername)
#    if ll[1]<200:
#        DICTlayers['chl'].append(layername)

Nlayers = {}
for var in varLIST:
    Nlayers[var] = len(DICTlayers[var])

DICTflagvar = {
    'chl': 0,
    'nit': 1,
}

Ndates = TLmis.nTimes
for misfile,datemis in zip(TLmis.filelist,TLmis.Timelist):
    date8 = datemis.strftime('%Y%m%d')
    #print misfile
    req = requestors.Daily_req(datemis.year,datemis.month,datemis.day)
    misALL = np.loadtxt(misfile,skiprows=1)
    LIST = {}
    for var in varLIST:
        LIST[var] = [i for i in range(1+2*Nlayers[var])]
        LIST[var][0] = datemis
        misvar = misALL[misALL[:,1]==DICTflagvar[var]]
    
        wmovar = set(misvar[:,-1])

        for il,ll in enumerate(DICTlayers[var]):
        # print '   ...' + ll
            LISTmeanmis = [[] for ii in range(2)] 
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
                    bias = np.nanmean(misint)
                    LISTmeanmis[0].append(rmsd)
                    LISTmeanmis[1].append(bias)
            if len(LISTmeanmis)>0:
                meanmis = np.nanmean(LISTmeanmis[0])
                meanbias = np.nanmean(LISTmeanmis[1])
            else:
                meanmis = np.nan
                meanbias = np.nan
            LIST[var][1+il] = meanmis
            LIST[var][1+Nlayers[var]+il] = meanbias


        nomefile = 'poststats.' + date8 + '.' + DICTvarname[var] + '.npy'
        np.save(OUTDIR + nomefile, LIST[var])







