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

    parser.add_argument('--maskfile', '-m',
                        type=str,
                        default=None,
                        required=True,
                        help="")

    return parser.parse_args()


args = argument()

import numpy as np
from commons.mask import Mask
from commons.submask import SubMask
from commons.dataextractor import DataExtractor
from commons.utils import addsep
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins import V2



INDADIR = addsep(args.indadir)
OUTDIR = addsep(args.outdir)
TheMask = Mask(args.maskfile)

_,jpj,jpi = TheMask.shape
mask200 = TheMask.mask_at_level(200)


Sup_mask = TheMask.cut_at_level(0)
dtype = [(sub.name, np.bool) for sub in V2.P]
npSub = {}
SUB = np.zeros((jpj,jpi),dtype=dtype)
Nsub = 1
for sub in V2.Pred:
    Nsub += 1
    sbmask         = SubMask(sub,maskobject=Sup_mask).mask
    SUB[sub.name]  = sbmask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]
    npSub[sub.name] = np.sum(SUB[sub.name])

npSub['med'] = np.sum(SUB['med'])


TLmis = TimeList.fromfilenames(None, INDADIR,  \
             '????????.chl_mis.nc', \
             prefix='',dateformat='%Y%m%d')


for misfile,datemis in zip(TLmis.filelist,TLmis.Timelist):
    date8 = datemis.strftime('%Y%m%d')
    misfmean = np.zeros((4,Nsub))
    misfmean[:] = np.nan
    De = DataExtractor(TheMask,filename=misfile,varname='misfchl',dimvar=2)
    chlmisf = De.values
    chlmisf[chlmisf>1.e+19] = np.nan


    for isub,sub in enumerate(V2.P):
        chlsub = chlmisf[SUB[sub.name]]
        if np.all(np.isnan(chlsub))==False:
            misfmean[0,isub] = (np.nanmean(chlsub**2))**.5
        misfmean[3,isub] = np.float(np.sum(np.isfinite(chlsub)))/npSub[sub.name]
        chlsub = chlmisf[SUB[sub.name] & (mask200==True)]
        if np.all(np.isnan(chlsub))==False:
            misfmean[1,isub] = (np.nanmean(chlsub**2))**.5
        chlsub = chlmisf[SUB[sub.name] & (mask200==False)]
        if np.all(np.isnan(chlsub))==False:
            misfmean[2,isub] = (np.nanmean(chlsub**2))**.5

    fileout = OUTDIR + date8 + '_meanmisf.npy'
    np.save(fileout,misfmean)


