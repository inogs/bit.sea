from commons.time_interval import TimeInterval

starttime='20160301'
end__time='20160308'
TI=TimeInterval(starttime,end__time,'%Y%m%d')

import basins.OGS as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
from commons.mask import Mask
from profiler import *

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')

Profilelist=bio_float.FloatSelector(None,TI,OGS.med)

M = Matchup_Manager(TI,INPUTDIR,BASEDIR)
M.getFloatMatchups(Profilelist,TheMask.zlevels,read_adjusted=False)

