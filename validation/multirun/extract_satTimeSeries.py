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
import pickle
import datetime

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
    for ii,infile in enumerate(filelist):
        datef = os.path.basename(infile)[0:8]
        print(datef)
        CH = NC.netcdf_file(infile,'r')
        chl = CH.variables[var].data.copy()
        for isub,sub in enumerate(SUB_LIST):
            print(sub)
            masksub = (SUB[sub]) & (chl>0.) & (mask200) #open sea
            chsub = chl[masksub]
            meanopen[isub,ii] = np.mean(chsub)
            stdopen[isub,ii] = np.std(chsub)
            masksub = (SUB[sub]) & (chl>0.) & (~mask200) #coast
            chsub = chl[masksub]
            meancoast[isub,ii] = np.mean(chsub)
            stdcoast[isub,ii] = np.std(chsub)
    stats_out[0,:] = meanopen
    stats_out[1,:] = meancoast
    stats_out[2,:] = stdopen
    stats_out[3,:] = stdcoast

    return stats_out;

fileslists = glob.glob(args.indir + args.year + '*')
fileslists.sort()
DATES = []
if fileslists:
   print('saving sat data')
   for ii in range(len(fileslists)):
       infile = fileslists[ii]
       datef = os.path.basename(infile)[0:8]
       T=datetime.datetime.strptime(datef,'%Y%m%d')
       DATES.append(T)

   stats_sat = stats_open_coast(fileslists,'lchlm')

   LISTsave = [i for i in range(3)]

   LISTsave[0] = DATES
   LISTsave[1] = stats_sat[0,:] # meansatopen
   LISTsave[2] = stats_sat[2,:] # stdsatopen

   outfile = args.outdir + 'satstats.open' + args.year + '.pkl'
   print(outfile)
   fid = open(outfile,'wb')
   pickle.dump(LISTsave,fid)
   fid.close()

   LISTsave = [i for i in range(3)]

   LISTsave[0] = DATES
   LISTsave[1] = stats_sat[1,:] # meansatcoast
   LISTsave[2] = stats_sat[3,:] # stdsatcoast

   outfile = args.outdir + 'satstats.coast' + args.year + '.pkl'
   print(outfile)
   fid = open(outfile,'wb')
   pickle.dump(LISTsave,fid)
   fid.close()


else: print('filelist sat empty!')
