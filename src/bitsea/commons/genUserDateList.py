from datetime import timedelta
from dateutil.relativedelta import relativedelta 
import datetime


def readTimeString(date17):

    year  = int(date17[0:4])
    month = int(date17[4:6])
    day   = int(date17[6:8])
    hr    = int(date17[9:11])
    mn    = int(date17[12:14])
    sec   = int(date17[15:  ])
    d = datetime.datetime(year, month, day,hr,mn,sec);
    return d;

def getCouplesTimeList(datestart,dateend,days=0, seconds=0, minutes=0, hours=0, years=0, months=0,overlap_hours=0):
    Datestart=readTimeString(datestart);
    Dateend  =readTimeString(dateend  );
    td = relativedelta(days=days,seconds=seconds, minutes=minutes, hours=hours, years=years,months=months)
    tdOverlap = relativedelta(hours=overlap_hours)
    tstart = Datestart
    tend__ = Datestart +td + tdOverlap
    TIMELIST=[]
        
    while tend__ <= Dateend:        
        TIMELIST.append([tstart,tend__])    
        tstart = tstart + td
        tend__ = tend__ + td #+tdOverlap
    if len(TIMELIST)>0:    
        lastcalculated=TIMELIST[-1][1]

        if lastcalculated < Dateend:
            TIMELIST.append([lastcalculated, Dateend ])
        
    return TIMELIST 

def getTimeList(datestart,dateend,days=0, seconds=0, minutes=0, hours=0, years=0, months=0):
    if isinstance(datestart,str):
        Datestart=readTimeString(datestart)
        Dateend  =readTimeString(dateend  )
    if isinstance(datestart, datetime.datetime):
        Datestart = datestart
        Dateend   = dateend   

    
    td = relativedelta(days=days,seconds=seconds, minutes=minutes, hours=hours, years=years,months=months)

    #print(td)
    tend  = Datestart+td
    TIMELIST=[Datestart]
        
    while tend <= Dateend:        
        TIMELIST.append(tend)            
        tend = tend + td

        
    return TIMELIST

if __name__ == '__main__':
    Days17 = getTimeList("19970127-00:00:00","19970601-00:00:00", days=17)
    Min15 = getTimeList("20180301-00:00:00","20200310-00:00:00", minutes=15)
    T = getCouplesTimeList("20190101-00:00:00", "20190501-00:00:00", months=1,  overlap_hours=12)










