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

def getCouplesTimeList(datestart,dateend,deltastring,overlapstring):
    Datestart=readTimeString(datestart);
    Dateend  =readTimeString(dateend  );
    td        = relativedelta(months = 1);
    tdOverlap = relativedelta(hours = 12) 
    tdcommand="td        =relativedelta(" +   deltastring + ")"; exec tdcommand    
    tdcommand="tdOverlap =relativedelta(" + overlapstring + ")"; exec tdcommand
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

def getTimeList(datestart,dateend,deltastring):
    if isinstance(datestart,str):
        Datestart=readTimeString(datestart)
        Dateend  =readTimeString(dateend  )
    if isinstance(datestart, datetime.datetime):
        Datestart = datestart
        Dateend   = dateend   
    td = relativedelta(months = 1); 
    tdcommand="td=relativedelta(" + deltastring + ")"
    exec tdcommand
    
    tend  = Datestart+td
    TIMELIST=[Datestart]
        
    while tend <= Dateend:        
        TIMELIST.append(tend)            
        tend = tend + td

        
    return TIMELIST











