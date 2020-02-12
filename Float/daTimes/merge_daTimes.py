import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    requires frequency of float DA
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--daTimes_sat',"-s",
                                type = str,
                                required = True,
                                help = 'daTimes_sat file')
    parser.add_argument(   '--daTimes_dir','-d',
                                type = str,
                                required = True,
                                help = 'output dir of profiles_dates_DAfreq.py')
    parser.add_argument(   '--outfile','-o',
                                type = str,
                                required = True,
                                help = 'daTimes used by ogstm.xx')
    return parser.parse_args()

args = argument()


from instruments.superfloat import FloatSelector
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import datetime
import numpy as np
from commons.utils import addsep

DIR=addsep(args.daTimes_dir)

LIStvarfloat = ['P_l']
LIStvarfloat = ['N3n','P_l']

# Merge with satellite DA dates 

dateDAsat = []
fid = open(args.daTimes_sat,'r')
for line in fid:
    dateDAsat.append(line.rstrip())

fid.close()

dateDAfloatvar = {}
for varfloat in LIStvarfloat:
    print varfloat
    dateDAfloatvar[varfloat] = []
    fid = open(DIR + 'daTimes_float_' + varfloat,'r')
    for line in fid:
        dateDAfloatvar[varfloat].append(line.rstrip())
    fid.close()

# removing dates for only N3n because ogsmt crashes there
dateDAfloatvar['N3n'] = [p for p in dateDAfloatvar['N3n'] if p in dateDAfloatvar['P_l']]

dateDAfloat = []
for varfloat in LIStvarfloat:
    dateDAfloat = list(set(dateDAfloat + dateDAfloatvar[varfloat]))

dateDAfloat.sort()
suffix = ''
for varfloat in LIStvarfloat:
    suffix = suffix + varfloat
np.savetxt(DIR + 'daTimes_float_' + suffix,dateDAfloat,fmt='%s')


dateDA = list(set(dateDAsat + dateDAfloat))

dateDA.sort()
np.savetxt(args.outfile,dateDA,fmt='%s')
