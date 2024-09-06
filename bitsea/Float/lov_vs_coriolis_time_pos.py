from instruments import lovbio_float
from instruments import bio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import numpy as np

TI = TimeInterval('2012','2020','%Y')
R = Rectangle(-6,36,30,46)
PROFILES_LOV=lovbio_float.FloatSelector(None, TI, R)

nLOVfiles=len(PROFILES_LOV)
BAD_PROFILES=[]
DIFF=np.zeros((nLOVfiles,3), np.float32)
for ip, p in enumerate(PROFILES_LOV):
    
    pCor = bio_float.from_lov_profile(p)
    
    if pCor is None:
        BAD_PROFILES.append(p)
    else:
        dt = p.time - pCor.time
        DIFF[ip,0]=p.lon - pCor.lon
        DIFF[ip,1]=p.lat - pCor.lat
        DIFF[ip,2]=dt.total_seconds()
        if ( (np.abs(DIFF[ip,0])>0.2) | (np.abs(DIFF[ip,1])>0.2) ):
            print "too much diff"
            

import matplotlib.pyplot as pl
pl.close('all')
fig,ax=pl.subplots()
ax.plot(DIFF[:,0], 'r.', label='lon')
ax.plot(DIFF[:,1], 'b.', label='lat')
ax.legend()
fig.show()

fig,ax=pl.subplots()
ax.plot(DIFF[:,2], 'g.', label='time')
ax.legend()
fig.show()
