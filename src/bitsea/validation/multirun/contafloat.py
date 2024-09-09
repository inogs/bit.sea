import os
from profiler_comparison2015 import *
from basins import OGS
import pylab as pl
#import matplotlib.dates as mdates
from timeseries.plot import *
from instruments.matchup_manager import Matchup_Manager
from instruments.var_conversions import LOVFLOATVARS
from instruments import lovbio_float as bio_float
from instruments import matchup_manager
from commons.utils import addsep
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import numpy as np


T_start = DATESTART
T_end = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start, '%Y%m%d'))
T_end2num = mpldates.date2num(datetime.strptime(T_end, '%Y%m%d'))


MM = Matchup_Manager(ALL_PROFILES, TL, BASEDIR['RSTaft'])


DICTadj = {
    'P_l': True,
    'N3n': True,
}


Profilelist_1 = bio_float.FloatSelector(None, TI1, OGS.med)
wmo_list = bio_float.get_wmo_list(Profilelist_1)

print wmo_list


timelabel_list = {}
for var in ['P_l','N3n']:
    timelabel_list[var] = []

for j in range(0, len(wmo_list)):
    print wmo_list[j]
#     if not(wmo_list[j]=='6902700'): continue
    list_float_track = bio_float.filter_by_wmo(Profilelist_1, wmo_list[j])
    print len(list_float_track)

    for ip, p in enumerate(list_float_track):
        timef12 = datetime(p.time.year,p.time.month,p.time.day,12,0)
        for varcheck in ['P_l','N3n']:
            _, Profile, _ = p.read(LOVFLOATVARS[varcheck],DICTadj[varcheck])
            if((Profile.size!=0) and (timef12 in TL.Timelist)):
                timelabel_list[varcheck].append(timef12)

#for var in ['P_l','N3n']:
#    timelabel_list[var].sort()


dates = {}
nprofiles_dates = {}
nprofiles_month = {}
TOTnp = {}
for var in ['P_l','N3n']:
    dates[var] = []
    nprofiles_dates[var] = []
    TOTnp[var] = 0
    nprofiles_month[var] = np.zeros(12)
    for dd in np.unique(timelabel_list[var]):
        dates[var].append(dd)
        nprof = 0
        for datep in timelabel_list[var]:
            if datep==dd: nprof+=1
        nprofiles_dates[var].append(nprof)
        TOTnp[var] += nprof
        indm = dd.month-1
        nprofiles_month[var][indm] += nprof


for var in ['P_l','N3n']: print TOTnp[var]

pl.close('all')

pl.figure(figsize=[9,5])
pl.bar(dates['P_l'],nprofiles_dates['P_l'],1, \
        label='P_l ' + np.str(TOTnp['P_l']))
pl.bar(dates['N3n'],nprofiles_dates['N3n'],1, \
        label='N3n ' + np.str(TOTnp['N3n']))

pl.grid()
pl.legend()

pl.savefig('/gpfs/scratch/userexternal/ateruzzi/' + \
           'ELAB_DAFloat/histNfloat.png')


pl.figure(figsize=[9,5])
pl.bar(range(1,13),nprofiles_month['P_l'],1, \
        label='P_l ' + np.str(TOTnp['P_l']))
pl.bar(range(1,13),nprofiles_month['N3n'],1, \
        label='N3n ' + np.str(TOTnp['N3n']))

pl.grid()
pl.legend()

pl.savefig('/gpfs/scratch/userexternal/ateruzzi/' + \
           'ELAB_DAFloat/histNfloat_month.png')

pl.show(block=False)



