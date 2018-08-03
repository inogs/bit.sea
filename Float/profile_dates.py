from instruments.lovbio_float import FloatSelector
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import datetime
import numpy as np


DATESTART = '20110101'
DATEEND = '20190112'

var = 'P_l' #'N3n' 'O2o'

TI = TimeInterval(DATESTART,DATEEND,'%Y%m%d')

ALL_PROFILES = FloatSelector(var,TI,Rectangle(-6,36,30,46))
# Verify if necessary to use a variable specification
# based on the 3dav namelist

datelist = []
for p in ALL_PROFILES:
    # print p.time
    tt = p.time.replace(hour=0,minute=0,second=0)
    datelist.append(datetime.datetime.strftime(tt,'%Y%m%d-%H:%M:%S'))


dateDAfloat = list(set(datelist))
dateDAfloat.sort()

np.savetxt('daTimes_float',dateDAfloat,fmt='%s')

# Merge with satellite DA dates 

dateDAsat = []
fid = open('daTimes_sat','r')
for line in fid:
    dateDAsat.append(line.rstrip())

fid.close()


dateDA = list(set(dateDAsat + dateDAfloat))

dateDA.sort()
