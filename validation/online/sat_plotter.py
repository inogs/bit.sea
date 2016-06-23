from sat_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
import pylab as pl

TI = TimeInterval('20160401','20160501','%Y%m%d')
ARCHIVE_DIR_WRONG="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2"
ARCHIVE_DIR      ="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2-dev"

A = timelistcontainer(TI,ARCHIVE_DIR_WRONG,'f0')
B = timelistcontainer(TI,ARCHIVE_DIR_WRONG,'f1')

times, v1= A.plotdata(A.bias, 'alb','open_sea')
times, v2= A.plotdata(B.number, 'alb','open_sea')

fig, ax = pl.subplots()
ax.plot(times,v2)
