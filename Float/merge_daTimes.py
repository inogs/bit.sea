from instruments.lovbio_float import FloatSelector
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import datetime
import numpy as np



# Merge with satellite DA dates 

dateDAsat = []
fid = open('daTimes_sat','r')
for line in fid:
    dateDAsat.append(line.rstrip())

fid.close()

dateDAfloat = []
fid = open('daTimes_float','r')
for line in fid:
    dateDAfloat.append(line.rstrip())

fid.close()


dateDA = list(set(dateDAsat + dateDAfloat))

dateDA.sort()
np.savetxt('daTimes',dateDA,fmt='%s')
