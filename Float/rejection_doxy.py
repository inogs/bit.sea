from basins.region import Rectangle
from commons.time_interval import TimeInterval
from instruments import bio_float
from Float import superfloat_generator
import commons.genUserDateList as DL
from commons.Timelist import TimeList
import numpy as np
import matplotlib.pyplot as pl
import sys

weekly=DL.getTimeList("20200107-00:00:00", "20210502-12:00:00", "days=7")
TL = TimeList(weekly)


var = 'DOXY'
R = Rectangle(-6,36,30,46)

#TI = TimeInterval('20210401','20210407','%Y%m%d')
WEEKLY_REQS=TL.getOwnList()


nWeeks = len(WEEKLY_REQS)

TOTAL = np.zeros((nWeeks), np.int)
REJ  = np.zeros((nWeeks), np.int)
for ireq, req in enumerate(WEEKLY_REQS):
    print ireq
    PROFILE_LIST=bio_float.FloatSelector(var, req.time_interval, R)
    rejected_list=[]
    for p in PROFILE_LIST:
        if p._my_float.status_var('DOXY') == 'R' :
            rejected_list.append(p)
        else:
            superfloat_name=p._my_float.filename.replace('CORIOLIS','SUPERFLOAT')
            if not superfloat_generator.exist_valid_variable('DOXY', superfloat_name):
                rejected_list.append(p)
    TOTAL[ireq] = len(PROFILE_LIST)
    REJ  [ireq] = len(rejected_list)

fig, ax = pl.subplots(figsize=(16,4))
ax.plot(weekly, TOTAL,'b', label='Coriolis')
ax.plot(weekly, REJ,  'r', label='rejected')
ax.plot(weekly, TOTAL-REJ,  'g', label='good')
ax.legend()
ax.grid()
fig.savefig('doxy_rejection.png')
