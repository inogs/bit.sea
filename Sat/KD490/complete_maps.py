import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    complete weekly maps using monthly mean or climatology at model resolution
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--maskfile', '-m',
                        type=str,
                        default=None,
                        required=True,
                        help=''' Path of maskfile''')

    parser.add_argument('--weeklydir', '-w',
                        type=str,
                        default=None,
                        required=True,
                        help=" Path of weekly files")

    parser.add_argument('--monthlydir', '-d',
                        type=str,
                        default=None,
                        required=True,
                        help=" Path of monthly files")

    parser.add_argument('--climafile', '-c',
                        type=str,
                        default=None,
                        required=True,
                        help=" File of filled climatology")

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help=" outdir")

    parser.add_argument('--outind', '-n',
                        type=str,
                        default=None,
                        required=True,
                        help=" outdir for indexes")

    return parser.parse_args()

args = argument()


import numpy as np
import glob,os
import datetime
import netCDF4 as NC4
import scipy.io.netcdf as NC
from commons.mask import Mask
from commons.utils import addsep
#from maskload import tmask,jpi,jpj

fmt = '%Y%m%d'


WEEKLY_DIR=addsep(args.weeklydir)
MONTHLY_DIR=addsep(args.monthlydir)
fileclim=args.climafile
OUTDIR = addsep(args.outdir)
OUTIND = addsep(args.outind)
TheMask = Mask(args.maskfile)

_,jpj,jpi = TheMask.shape

tmaskS = TheMask.mask[0,:,:]

filelist = glob.glob(WEEKLY_DIR + '/*_d-OC_CNR-L3-KD490-MedOC4AD4_SAM_1KM-MED-REP*.nc')
filelist.sort()

Kc = NC4.Dataset(fileclim,'r') 

# indexes to recontruct maps
nmask = np.sum(tmaskS)
indI = np.zeros(nmask,dtype=int)
indJ = np.zeros(nmask,dtype=int)
ip = 0
for jj in range(jpj):
    for ii in range(jpi):
        if tmaskS[jj,ii]:
           indI[ip] = ii
           indJ[ip] = jj
           ip = ip+1
npfile = OUTIND + '/indI.npy'
np.save(npfile,indI)
npfile = OUTIND + '/indJ.npy'
np.save(npfile,indJ)

iclim = 0
for infile in filelist:
    datec = os.path.basename(infile)[0:8]
    print(datec)
    K = NC4.Dataset(infile,'r')
    #kext = np.array(K.variables['KD490'][0,:,:])
    kext = np.array(K.variables['KD490'][:,:])
    maskk = (kext == -999.) & (tmaskS==1)
    nump = np.sum(maskk)
    print(nump)
    mstr = os.path.basename(infile)[0:6]
    filemonth = MONTHLY_DIR + '/' + mstr + '_d-OC_CNR-L3-KD490-MedOC4AD4_SAM_1KM-MED-REP-v01.nc'
    Km = NC4.Dataset(filemonth,'r')
    # kmnth = np.array(Km.variables['KD490'][0,:,:])
    kmnth = np.array(Km.variables['KD490'][:,:])
    kext[maskk] = kmnth[maskk]

    mask_cntrl = (kext == -999.) & (tmaskS==1)
    ncntrl = np.sum(mask_cntrl)
    if ncntrl>0 :
       print('need climatology')
       iclim = iclim+1
    #    fileclim = CLIMA_DIR + '/KD490_Climatology_24_filled.nc'
       fnex = 1
       dt = datetime.datetime.strptime(datec, fmt)
       tt = dt.timetuple()
       julian_day=tt.tm_yday
       if julian_day==366: julian_day=365
       kclim = np.array(Kc.variables['Mean'][julian_day-1,:,:])
       kclim[kclim>1.e+19] = -999.
       kext[mask_cntrl] = kclim[mask_cntrl]

    kext[tmaskS==0] = -999.
    filenp = OUTDIR +'/' + os.path.basename(infile)[0:8] + '_compl.npy'
    #filenp = OUTDIR +'DAYS7_16compl/' + os.path.basename(infile)[0:8] + '_compl.npy'
    np.save(filenp,kext)
