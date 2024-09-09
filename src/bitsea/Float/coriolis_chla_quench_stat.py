from instruments import bio_float
#from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import matplotlib.pyplot as pl
import sys
from commons import calculated_depths
import numpy as np
import scipy.io.netcdf as NC
TI = TimeInterval('2012','2020','%Y')
R = Rectangle(-6,36,30,46)
PROFILES_COR =bio_float.FloatSelector('CHLA', TI, R)

wmo_list = bio_float.get_wmo_list(PROFILES_COR)


LINES=[]
fmt="%s\t"*6+ "\n"
line = fmt %( 'flag', 'month', 'MLD', 'Pres_top', 'quench_value', 'ID')
LINES.append(line)

for iwmo, wmo in enumerate(wmo_list[57:]):
    print iwmo, wmo
    if wmo == "6901765" : continue # save time
    track_list = bio_float.filter_by_wmo(PROFILES_COR, wmo)
    for pCor in track_list:
        if pCor._my_float.status_var('CHLA')=='A':
            PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=True)
            if len(ValueC) < 5:
                print "few values in Coriolis for " + pCor._my_float.filename
                continue
            PresT, Temp, _ = pCor.read('TEMP', read_adjusted=False)
            mld = calculated_depths.mld(Temp, PresT, zref=0, deltaTemp=0.1)
            
            quench_value = np.nan
            if mld < 80:
                if PresC[0]<= mld+5:
                    ii = PresC<= mld + 5
                    if ii.sum()==0:
                        print "found zero"
                        sys.exit()
                    descending_ordered = np.sort(ValueC[ii])[-1::-1]
                    n_chosen = min(3, len(descending_ordered))
                    if n_chosen==0: n_chosen = 1
                    quench_value = descending_ordered[:n_chosen].mean()
                    flag=3
                else:
                    flag=2
            else:
                flag=1
            fmt = "%d\t"*2 + "%4.2f\t"*3 + "%s\n"
            line=fmt  %(flag, pCor.time.month, mld, PresC[0], quench_value, pCor.ID())
            LINES.append(line)
fid = open('Quench_statistics.txt','w')
fid.writelines(LINES)
fid.close()


dtype=[('flag', np.int), ('month', np.int), ('MLD', np.float32), ('Pres_top', np.float32), ('quench_value', np.float32), ('ID', "S200")]