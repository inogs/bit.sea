import argparse
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

    parser.add_argument(   '--hourda',"-t",
                                type = str,
                                required = True,
                                help = 'time during the day at which DA should be performed')
    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--outfile','-o',
                                type = str,
                                required = True,
                                help = 'path of the output text file ')
    return parser.parse_args()

args = argument()

from commons.time_interval import TimeInterval
from commons import timerequestors
from instruments import superfloat as biofloat
from instruments.var_conversions import FLOATVARS
import argparse
import basins.OGS as OGS
import commons.genUserDateList as DL
import datetime
import numpy as np

DAfreq = np.int(args.dafreq) # days
varMODEL = args.varda
hourDA = np.int(args.hourda)

deltatimeDA = datetime.timedelta(DAfreq)

DATESTART = args.datestart
DATEEND = args.dateend


var = FLOATVARS[varMODEL] #'N3n' 'O2o'

read_adjusted = {
    FLOATVARS['P_l']: True,
    FLOATVARS['N3n']: True,
}


TL = DL.getTimeList(DATESTART + '-00:00:00',DATEEND + '-00:00:00', \
                    'days=' + np.str(DAfreq))


DAfloatdates = []
NnoDAdates = 0
for dateref in TL:
    datefreq = timerequestors.Interval_req(dateref.year,dateref.month,dateref.day, \
                        days=DAfreq)
    dateend = datefreq.time_interval.end_time - datetime.timedelta(days=1)
    datefreq.time_interval.end_time = datetime.datetime(dateend.year, \
                                                        dateend.month, \
                                                        dateend.day, 23, 59)
    PROFILESdateref = biofloat.FloatSelector(var,datefreq.time_interval,OGS.med)

    Goodlist = []
    WMOlist = biofloat.get_wmo_list(PROFILESdateref)
    for wmo in WMOlist:
        SubProfilelist_1 = biofloat.filter_by_wmo(PROFILESdateref,wmo)
        for i in SubProfilelist_1:
            _, Profile, _ = i.read(var,read_adjusted[var])   #Profile.shape,Profile.size, np.mean(Profile)
            if(Profile.size!=0) : Goodlist.append(i)

        if (Goodlist!=[]):
            # print(wmo)
            break

    if (Goodlist!=[]):
        dateDA = datetime.datetime(dateref.year, \
                dateref.month,dateref.day,hourDA,00)
        print dateDA
                                   
        DAfloatdates.append(datetime.datetime.strftime(dateDA,'%Y%m%d-%H:%M:%S'))
    else:
        NnoDAdates += 1


np.savetxt(args.outfile,DAfloatdates,fmt='%s')

