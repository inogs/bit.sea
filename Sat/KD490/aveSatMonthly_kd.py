import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Fill climatology with nearest
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--checkdir', '-i',
                                type = str,
                                required =True,
                                help = ''' Climatology input .nc'''
                                )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'Output dir ')
    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')



    return parser.parse_args()

args = argument()

from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.utils import addsep
from postproc import masks
import numpy as np
import os
import Sat.SatManager as Sat

#CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/CHECKED/"
CHECKDIR = addsep(args.checkdir)
#MONTHLYDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/MONTHLY/ORIGMESH/"
MONTHLYDIR = addsep(args.outdir)

maskSat = getattr(masks,args.mesh)

reset = True

Timestart="19990101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
#IonamesFile = '../../postproc/IOnames_sat.xml'
#IOname = IOnames.IOnames(IonamesFile)
suffix = os.path.basename(TLCheck.filelist[0])[8:]
MONTHLY_reqs=TLCheck.getMonthlist()

jpi = maskSat.jpi
jpj = maskSat.jpj


for req in MONTHLY_reqs:
    outfile = req.string + suffix
    outpathfile = MONTHLYDIR + outfile
    conditionToSkip = (os.path.exists(outpathfile)) and (not reset)

    if conditionToSkip: continue

    print outfile
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        print '  ... %s of %s' %(iFrame+1,len(ii))
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile,'KD490')
        M[iFrame,:,:] = CHL
    print ' average ...'
    Kext_OUT = Sat.averager(M)
    print ' write ' + outpathfile
    Sat.dumpGenericNativefile(outpathfile, Kext_OUT,"KD490",mesh=maskSat)
