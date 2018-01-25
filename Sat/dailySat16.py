from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons import IOnames
import numpy as np
import argparse
import os
import SatManager as Sat

def argument():
    parser = argparse.ArgumentParser(description = '''
    Apply check based on climatology to sat ORIG files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--checkdir', '-i',
                                type = str,
                                required = True,
                                help = ''' CHECKED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/'''

                                )

    parser.add_argument(   '--weeklydir', '-o',
                                type = str,
                                required = True,
                                help = ''' CHECKED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/WEEKLY/'''

                                )


    return parser.parse_args()

args = argument()



CHECKDIR = args.checkdir
print('anna')
print(CHECKDIR)
WEEKLYDIR= args.weeklydir
print('anna')
print(WEEKLYDIR)

reset = False

Timestart="19500101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
IonamesFile = '../postproc/IOnames_sat.xml'
IOname = IOnames.IOnames(IonamesFile)


jpi = Sat.NativeMesh.jpi
jpj = Sat.NativeMesh.jpj


for infile in TLCheck.filelist:
    print(infile)
    print(os.path.basename(infile))
    outfile = os.path.basename(infile)
    outpathfile = WEEKLYDIR + outfile
    conditionToSkip = (os.path.exists(outpathfile)) and (not reset)

    if conditionToSkip: continue

    M = np.zeros((jpj,jpi),np.float32)
    CHL = Sat.readfromfile(infile)
    M[:,:] = CHL
    CHL_OUT = M.copy()
    Sat.dumpV4file(outpathfile, CHL_OUT, varname='CHL')

