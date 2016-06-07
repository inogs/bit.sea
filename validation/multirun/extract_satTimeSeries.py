import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Extracts statistcs for sat
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required =True,
                                help = ''' Output dir'''
                                )

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                required =True,
                                help = ''' Input dir'''
                                )

    parser.add_argument(   '--year', '-y',
                                type = str,
                                required =True,
                                help = ''' Year of evaluation'''
                                )

    return parser.parse_args()

args = argument()


import numpy as np
import scipy.io.netcdf as NC
import glob,os

from commons.mask import Mask
from commons.submask import SubMask
from basins import OGS



maskfile = os.getenv('MASKFILE')
TheMask = Mask(maskfile)
tk_m     = TheMask.getDepthIndex(200)
mask200   = TheMask.mask_at_level(200)

jpk,jpj,jpi =TheMask.shape
dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
for sub in OGS.P:
    sbmask         = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name]  = sbmask[0,:,:]

SUB_LIST = ['alb', 'sww', 'swe', 'nwm', 'tyr', 'adn', 'ads', 'aeg', 'ion', 'lev', 'med']
numsub = len(SUB_LIST)

def stats_open_coast(filelist,var):
    numd = len(filelist)
    print('Reading... ' + filelist[0] + ' -> ' + filelist[-1])
    stats_out = np.zeros((4,numsub,numd))
    meanopen = np.zeros((numsub,numd))
    meancoast = np.zeros((numsub,numd))
    stdopen = np.zeros((numsub,numd))
    stdcoast = np.zeros((numsub,numd))
    DATES = []
    for ii in range(numd):
        infile = filelist[ii]
        datef = os.path.basename(infile)[5:12]
        DATES.append(datef)
        CH = NC.netcdf_file(infile,'r')
        chl = CH.variables[var].data.copy()
        for isub in range(len(SUB_LIST)):
            sub = SUB_LIST[isub]
            masksub = (SUB[sub][0,:]) & (chl<1.e+19) & (mask200) #open sea
            chsub = chl[masksub]
            meanopen[isub,ii] = np.mean(chsub)
            stdopen[isub,ii] = np.std(chsub)
            masksub = (SUB[sub][0,:]) & (chl<1.e+19) & (~mask200) #coast
            chsub = chl[masksub]
            meancoast[isub,ii] = np.mean(chsub)
            stdcoast[isub,ii] = np.std(chsub)
    stats_out[0,:] = meanopen
    stats_out[1,:] = meancoast
    stats_out[2,:] = stdopen
    stats_out[3,:] = stdcoast
    return stats_out


fileslists = glob.glob(args.indir + args.year + '*')
fileslists.sort()
DATES = []
if fileslists:
   for ii in range(len(fileslists)):
       infile = fileslists[ii]
       datef = os.path.basename(infile)[0:8]
       DATES.append(datef)

   stats_sat = stats_open_coast(fileslists,'lchlm')
   filedat = args.outdir + 'datesat.npy'
   np.save(filedat,DATES)

npyfile = args.outdir + 'meansatopen' + args.year + '.npy'
np.save(npyfile,stats_sat[0,:])
npyfile = args.outdir + 'meansatcoast' + args.year + '.npy'
np.save(npyfile,stats_sat[1,:])
npyfile = args.outdir + 'stdsatopen' + args.year + '.npy'
np.save(npyfile,stats_sat[2,:])
npyfile = args.outdir + 'stdsatcoast' + args.year + '.npy'
np.save(npyfile,stats_sat[3,:])

