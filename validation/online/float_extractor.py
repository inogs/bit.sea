from commons.time_interval import TimeInterval

starttime='20160301'
end__time='20160308'
TI=TimeInterval(starttime,end__time,'%Y%m%d')

import basins.OGS as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
#from commons.layer import Layer

modelvarname='O2o'

Profilelist=bio_float.FloatSelector(None,TI,OGS.med)

WMO=set()
for p in Profilelist:
    F=p._my_float
    WMO.add(F.wmo)


#M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)


