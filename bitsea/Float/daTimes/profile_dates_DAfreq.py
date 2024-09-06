import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    requires frequency of float DA
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--dafreq',"-f",
                                type = str,
                                required = True,
                                help = 'frequency for float DA in days')
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

from commons import timerequestors
from instruments.var_conversions import FLOATVARS
import basins.OGS as OGS
import commons.genUserDateList as DL
import datetime
import numpy as np
import os,sys

profilesource=os.getenv("PROFILES_SOURCE")
if profilesource is None:
    print("Error: Environment variables PROFILES_SOURCE - superfloat or ppcon - must be defined.")
    sys.exit(1)
assert profilesource in ["superfloat", "ppcon"]
if profilesource=="superfloat":
    from instruments import superfloat as biofloat
if profilesource=="ppcon":
    from instruments import float_ppcon as biofloat



DAfreq = int(args.dafreq)
varMODEL = args.varda
hourDA = int(args.hourda)

deltatimeDA = datetime.timedelta(DAfreq)

DATESTART = args.datestart
DATEEND = args.dateend


var = FLOATVARS[varMODEL]
TL = DL.getTimeList(DATESTART + '-00:00:00',DATEEND + '-00:00:00', days=DAfreq)

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
            _, Profile, _ = i.read(var)
            if(Profile.size!=0) : Goodlist.append(i)

        if (Goodlist!=[]):
            break

    if (Goodlist!=[]):
        dateDA = datetime.datetime(dateref.year, \
                dateref.month,dateref.day,hourDA,00)
        print (dateDA)
                                   
        DAfloatdates.append(datetime.datetime.strftime(dateDA,'%Y%m%d-%H:%M:%S'))
    else:
        NnoDAdates += 1


np.savetxt(args.outfile,DAfloatdates,fmt='%s')
