from instruments.lovbio_float import FloatSelector
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import datetime
import numpy as np


LIStvarfloat = ['P_l']
LIStvarfloat = ['N3n','P_l']

# Merge with satellite DA dates 

dateDAsat = []
fid = open('daTimes_sat','r')
for line in fid:
    dateDAsat.append(line.rstrip())

fid.close()

dateDAfloatvar = {}
for varfloat in LIStvarfloat:
    print varfloat
    dateDAfloatvar[varfloat] = []
    fid = open('daTimes_float_' + varfloat,'r')
    for line in fid:
        dateDAfloatvar[varfloat].append(line.rstrip())
    fid.close()

dateDAfloat = []
for varfloat in LIStvarfloat:
    dateDAfloat = list(set(dateDAfloat + dateDAfloatvar[varfloat]))

dateDAfloat.sort()
suffix = ''
for varfloat in LIStvarfloat:
    suffix = suffix + varfloat
np.savetxt('daTimes_float_' + suffix,dateDAfloat,fmt='%s')


dateDA = list(set(dateDAsat + dateDAfloat))

dateDA.sort()
np.savetxt('daTimes',dateDA,fmt='%s')
