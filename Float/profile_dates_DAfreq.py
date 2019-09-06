from commons.time_interval import TimeInterval
from commons import timerequestors
from instruments.lovbio_float import FloatSelector
from instruments.lovbio_float import get_wmo_list
from instruments.lovbio_float import filter_by_wmo
from instruments.var_conversions import LOVFLOATVARS
import argparse
import basins.OGS as OGS
import commons.genUserDateList as DL
import datetime
import numpy as np

def argument():
    parser = argparse.ArgumentParser(description = '''
    requires frequency of float DA
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--dafreq',"-f",
                                type = str,
                                required = True,
                                help = 'frequency for float DA')

    parser.add_argument(   '--varda',"-v",
                                type = str,
                                required = True,
                                help = 'variable (model name) for float DA')

    return parser.parse_args()

args = argument()



DAfreq = np.int(args.dafreq) # days
varMODEL = args.varda

deltatimeDA = datetime.timedelta(DAfreq)

DATESTART = '20150102'
DATEEND = '20160101'


var = LOVFLOATVARS[varMODEL] #'N3n' 'O2o'

read_adjusted = {
    LOVFLOATVARS['P_l']: True,
    LOVFLOATVARS['N3n']: True,
}


TL = DL.getTimeList(DATESTART + '-00:00:00',DATEEND + '-00:00:00', \
                    'days=' + np.str(DAfreq))


DAfloatdates = []
NnoDAdates = 0
for dateref in TL[1:]:
    datefreq = timerequestors.Interval_req(dateref.year,dateref.month,dateref.day, \
                        days=DAfreq)
    dateend = datefreq.time_interval.end_time - datetime.timedelta(days=1)
    datefreq.time_interval.end_time = datetime.datetime(dateend.year, \
                                                        dateend.month, \
                                                        dateend.day, 23, 59)
    PROFILESdateref = FloatSelector(var,datefreq.time_interval,OGS.med)

    Goodlist = []
    WMOlist = get_wmo_list(PROFILESdateref)
    for wmo in WMOlist:
        SubProfilelist_1 = filter_by_wmo(PROFILESdateref,wmo)
        for i in SubProfilelist_1:
            _, Profile, _ = i.read(var,read_adjusted[var])   #Profile.shape,Profile.size, np.mean(Profile)
            if(Profile.size!=0) : Goodlist.append(i)

        if (Goodlist!=[]):
            # print(wmo)
            break
    
    if (Goodlist!=[]):
        print dateref
        DAfloatdates.append(datetime.datetime.strftime(dateref,'%Y%m%d-%H:%M:%S'))
    else:
        NnoDAdates += 1


np.savetxt('daTimes_floatfreq' + np.str(DAfreq) + '_' + varMODEL ,DAfloatdates,fmt='%s')

