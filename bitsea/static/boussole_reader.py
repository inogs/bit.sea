from commons import netcdf4
from datetime import datetime, timedelta
from commons.time_interval import TimeInterval
import numpy as np

times_seconds=netcdf4.readfile('out.nc', 'time')
depth= netcdf4.readfile('out.nc', 'depth')
chla = netcdf4.readfile('out.nc', 'chla')
lam  = netcdf4.readfile('out.nc', 'lambda')

Dref = datetime(1970,1,1,0,0,0)

unique_times, inverse, counts = np.unique(times_seconds,return_inverse=True, return_counts=True )

TIMELIST= [Dref + timedelta(seconds=s) for s in unique_times ]


TI = TimeInterval('200302','200304','%Y%m')

class my_profile():
    def __init__(self,time,depth,chla,lambda_matrix):
        self.time = time
        self.depth = depth
        self.chla = chla
        self.lambda_matrix= lambda_matrix

Profilelist=[]
for it,t in enumerate(TIMELIST):
    if TI.contains(t):
        ii=inverse==it
        p = my_profile(t,depth[ii],chla[ii],lam[ii,:])
        Profilelist.append(p)
