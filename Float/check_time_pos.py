from instruments import lovbio_float
from basins.region import Rectangle
from commons.time_interval import TimeInterval
import datetime
    
TI = TimeInterval('20100101','2050101','%Y%m%d')
R = Rectangle(-6,36,30,46)

PROFILE_LIST=lovbio_float.FloatSelector(None, TI, R)

WMOS=lovbio_float.get_wmo_list(PROFILE_LIST)

for wmo in WMOS:
    
    Profile_list = lovbio_float.filter_by_wmo(PROFILE_LIST, wmo)
    #first=Profile_list[0]
    #d=datetime.datetime(firsrt.time)
    
    TIMELIST=[p.time for p in Profile_list]
    nProfiles = len(Profile_list)
    for ip in range(1,nProfiles):
        if TIMELIST[ip] <= TIMELIST[ip-1]:
            print "PROBLEM in ", wmo, ip
            for k in range(-1,2):
                IP = ip+k
                print IP, Profile_list[IP]._my_float.filename, TIMELIST[IP]

    
