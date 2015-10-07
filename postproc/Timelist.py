import requestors
import genUserDateList as DL
import os,glob
import datetime
import numpy as np
import IOnames

def isInWindow( s1, e1, s2, e2):
    overlap = min(e1,e2) - max(s1,s2)
    res     = np.double( overlap.days*86400 + overlap.seconds )
    return res > 0
def overlapTime(s1,e1,s2,e2):

#------------------s1++++++++++e1--------------------- # Window of temporal average of model output file
#                              | 
#--------------------------s2++++++++++e2------------- # Window of aggregation defined in PLOT_LIST
#                           |  |
#--------------------------<++++>--------------------- # overlapping time window = theWindow

    theWindow = min(e1,e2) - max(s1,s2)     

# Convert time window in seconds and nornamlize over the aggreation time window
    res       = np.double( theWindow.days*86400 + theWindow.seconds ) 

# return only non-negative results

    return max( 0, res)

def computeTimeWindow(freqString,currentDate):

    if (freqString == 'daily'):   req = requestors.Daily_req(currentDate.year,currentDate.month,currentDate.day)
    if (freqString == 'weekly'):  req = requestors.Weekly_req(currentDate.year, currentDate.month,currentDate.day)
    if (freqString == 'monthly'): req = requestors.Monthly_req(currentDate.year, currentDate.month)

    return req.starttime, req.endtime

def getSeason(datetime_obj):
    '''
    Returns an integer indicating the season
    
    
    Assumption for integers indicating seasons:
    Winter = 0
    Spring = 1
    Summer = 2
    Fall   = 3
    '''
    if datetime_obj.month in [1,2,3] : return 0
    if datetime_obj.month in [4,5,6] : return 1
    if datetime_obj.month in [7,8,9] : return 2
    return 3

class TimeList():
    '''
    HYPOTHESIS: 
     - the Input directory has files containing a date in their names
     - every file refers to a single time frame
     - The date of in the file name can be considered as centered in their period
    The generated datetime list has all the files concerning the period indicated by datestart, dateend. 
    Some of these files can have the centered date out of that period.
    
    Example:
    
    INPUTDIR="/pico/scratch/userexternal/gbolzon0/Carbonatic/wrkdir/MODEL/AVE_FREQ_1/"
    TL = TimeList('20141201-00:00:00','20150701-00:00:00', INPUTDIR,"ave*N1p.nc", 'IOnames.xml')
    '''       
    def __init__(self, datestart,dateend, inputdir,searchstring,IOnamesfile):
     
        IOname = IOnames.IOnames(IOnamesfile)
        self.datestart     = datetime.datetime.strptime(datestart,IOname.Input.dateformat)
        self.dateend       = datetime.datetime.strptime(dateend  ,IOname.Input.dateformat)
        self.inputdir      = inputdir
        self.searchstring  = searchstring
        
        filelist_ALL = glob.glob(self.inputdir + self.searchstring)
        self.filelist=[]
        self.Timelist=[]
        External_filelist=[]
        External_timelist=[]
        for pathfile in filelist_ALL:
            filename   = os.path.basename(pathfile)
            datestr     = filename[IOname.Input.date_startpos:IOname.Input.date_endpos]
            actualtime = datetime.datetime.strptime(datestr,IOname.Input.dateformat)
            if actualtime >= self.datestart and actualtime <= self.dateend:        
                self.filelist.append(pathfile)
                self.Timelist.append(actualtime)
            else:
                External_filelist.append(pathfile)
                External_timelist.append(actualtime)
        
        self.filelist.sort()
        self.Timelist.sort()        
        self.inputFrequency= self.__searchFrequency__()
        
        if self.inputFrequency is not None:
            for iFrame, t in enumerate(External_timelist):
                s1,e1 = computeTimeWindow(self.inputFrequency, t)
                if isInWindow(s1, e1, self.datestart, self.dateend):
                    self.filelist.append(External_filelist[iFrame])
                    self.Timelist.append(External_timelist[iFrame])
            self.filelist.sort()
            self.Timelist.sort()
        self.nTimes   =len(self.filelist)
        
        if self.inputFrequency == 'daily':
            for iFrame, t in enumerate(self.Timelist):
                if t.hour==0:
                    newt = datetime.datetime(t.year,t.month,t.day,12,0,0)
                    self.Timelist[iFrame] = newt


    def __searchFrequency__(self):
        if len(self.Timelist)<2:
            timestr = self.datestart.strftime(" between %Y%m%d and ") +  self.dateend.strftime("%Y%m%d")
            print "Frequency cannot be calculated in " + self.inputdir + timestr
            return None
        mydiff = self.Timelist[1]-self.Timelist[0]
        if mydiff.days == 1:
            return "daily"
        if (mydiff.days) > 6 & (mydiff.days < 8 ) : 
            return "weekly" #centered in " + str(self.Timelist[1].isoweekday())
        if (mydiff.days > 28) & (mydiff.days < 32) :
            return "monthly"
        
    
    def __generaltimeselector__(self,requestor):
            SELECTION=[]
            weights  =[]
            
            if self.inputFrequency == "daily":
                for t in self.Timelist:
                    if (t>=requestor.starttime)  & (t<=requestor.endtime):
                        SELECTION.append(t)
                        weights.append(1.)
            if self.inputFrequency in ['weekly','monthly']:
                for t in self.Timelist:
                    s1,e1 = computeTimeWindow(self.inputFrequency,t);
                    s2    = requestor.starttime
                    e2    = requestor.endtime
                    weight = overlapTime(s1,e1,s2,e2)
                    if (weight > 0. ) : 
                        SELECTION.append(t)
                        weights.append(weight)
            return SELECTION , np.array(weights)           
        
    def select(self,requestor):
        '''
        indexes, weights = select(requestor)
        Returned values: 
         - a list of indexes (integers) indicating to access selected times (or files)
         - a numpy array of weights
        
         
        '''
            
        if isinstance(requestor, requestors.Monthly_req):
            SELECTION=[]
            weights = []
                           
            if self.inputFrequency == 'daily':
                for it, t in enumerate(self.Timelist):
                    if (t.year==requestor.year) & (t.month==requestor.month):
                        SELECTION.append(it)
                        weights.append(1.)
            if self.inputFrequency == 'weekly':
                for it,t in enumerate(self.Timelist):
                    s1,e1 = computeTimeWindow("weekly",t);
                    s2    = requestor.starttime
                    e2    = requestor.endtime
                    weight = overlapTime(s1,e1,s2,e2); 
                    if (weight > 0. ) : 
                        SELECTION.append(it)
                        weights.append(weight)
                    
            return SELECTION , np.array(weights)

        
        if isinstance(requestor, requestors.Weekly_req):
            assert self.inputFrequency != "monthly"
            assert self.inputFrequency != "weekly"
            
            SELECTION=[]
            weights  =[]
            
            if self.inputFrequency == "daily":
                for it,t in enumerate(self.Timelist):
                    if (t>=requestor.starttime) & (t<=requestor.endtime):
                        SELECTION.append(it)
                        weights.append(1.)                     
            return SELECTION , np.array(weights)
        
        
        if isinstance(requestor, requestors.Season_req):
            return self.__generaltimeselector__(requestor)
                
        if isinstance(requestor, requestors.Yearly_req):
            return self.__generaltimeselector__(requestor)
        
        if isinstance(requestor, requestors.Decadal_req):
            return self.__generaltimeselector__(requestor)
        
        if isinstance(requestor, requestors.Clim_day):
            assert self.inputFrequency == 'daily'
            SELECTION=[]
            weights  =[]
            for it,t in enumerate(self.Timelist):
                if (t.month == requestor.month) & (t.day == requestor.day):
                    SELECTION.append(it)
                    weights.append(1.)
            return SELECTION , np.array(weights)

        
        if isinstance(requestor,requestors.Clim_month):
            SELECTION = []
            weights   = []
            YEARLIST=self.getYearlist()
            for year_req in YEARLIST:
                req = requestors.Monthly_req(year_req.year, requestor.month)
                s,w = self.select(req)
                SELECTION.extend(s)
                weights.extend(w)
            
            return SELECTION , np.array(weights)
            
        if isinstance(requestor,requestors.Clim_season):
            SELECTION = []
            weights   = []
            YEARLIST=self.getYearlist()
            for year_req in YEARLIST:
                req = requestors.Season_req(year_req.year, requestor.season)
                s,w = self.select(req)
                SELECTION.extend(s)
                weights.extend(w)
            return SELECTION , np.array(weights)           
            
            
    def getWeeklyList(self,weekday):
        '''
        
        WeekList = getWeekList(weekday)
        Weekday is the same of the datetime.isoweekay() method.
        Monday == 1 ... Sunday == 7
        
        Returns an ordered list of requestors, Weekly_req objects.
        '''

        PossibleShifts=[3,2,1,-0,-1,-2,-3]
        for counter,day in enumerate(range(weekday-3,weekday+4)):
            interested_weekday = (day)%7
            if interested_weekday ==0 : interested_weekday = 7
            if interested_weekday == self.datestart.isoweekday():
                index = counter

        starting_centered_day = self.datestart + datetime.timedelta(days=PossibleShifts[index])

        TL=DL.getTimeList(starting_centered_day,self.Timelist[-1] , "days=7")
        REQ_LIST=[]
        for t in TL:
            m = requestors.Weekly_req(t.year,t.month,t.day)
            indexes,_ = self.select(m)
            if len(indexes)>0:
                REQ_LIST.append(m)
        return REQ_LIST             
            
    def getMonthlist(self, extrap=False):
        '''Returns an ordered list of requestors, Monthly_req objects. 
       By setting extrap=True, this method extrapolates out of the indicated period (Starttime,EndTime)  
       Example: if the input is weekly, centered in 20120301, and Starttime=20120301, 
       then extrapolation will return also 201202, because of the part of the week before the centered time. 
        ''' 
        t=self.Timelist[0]
        MONTH_LIST=[(t.year,t.month)]
        for t in self.Timelist:
            newmonth = (t.year, t.month)
            if newmonth not in MONTH_LIST: MONTH_LIST.append(newmonth)
        if extrap:
            pass
        else:
            MONTH_LIST_RED=[]
            firstMonth = datetime.datetime(self.datestart.year,self.datestart.month,1)
            lastMonth  = datetime.datetime(self.dateend.year,  self.dateend.month  ,1)
            for month in MONTH_LIST:
                firstOfMonth = datetime.datetime(month[0],month[1],1)
                if (firstOfMonth >= firstMonth) & (firstOfMonth <=lastMonth):
                    MONTH_LIST_RED.append(month) 
            MONTH_LIST = MONTH_LIST_RED
                
                
        
        REQ_LIST=[]
        for month in MONTH_LIST:
            m = requestors.Monthly_req(month[0],month[1])
            REQ_LIST.append(m)
        
        return REQ_LIST
    
    
    def getYearlist(self):
        '''
        Returns an list of requestors, Yearly_req objects.
        '''        
        YEARLIST=[self.datestart.year]
        for t in self.Timelist:
            if t.year not in YEARLIST: YEARLIST.append(t.year)
        REQ_LIST=[]
        for year in YEARLIST:
            m = requestors.Yearly_req(year)
            REQ_LIST.append(m)
        return REQ_LIST
    

    
    def getSeasonList(self,extrap=False):
        '''Returns an ordered list of requestors, Season_req objects. 
       By setting extrap=True, this method extrapolates out of the indicated period (Starttime,EndTime)  
       Example: if the input is weekly, centered in 20120301, and Starttime=20120301, 
       then extrapolation will return also 201202, because of the part of the week before the centered time. 
        '''         
        t=self.Timelist[0]
        SEASON_LIST=[(t.year, getSeason(t))]
        for t in self.Timelist:
            newSeason = (t.year, getSeason(t))
            if newSeason not in SEASON_LIST: SEASON_LIST.append(newSeason)
        if extrap:
            pass
        else:
            SEASON_LIST_RED=[]
            firstSeason=requestors.Season_req(self.datestart.year, getSeason(self.datestart))
            lastSeason = requestors.Season_req(self.dateend.year,  getSeason(self.dateend))
            for season in SEASON_LIST:
                req = requestors.Season_req(season[0],season[1])
                if (req.starttime >= firstSeason.starttime) & (req.endtime <=lastSeason.endtime):
                    SEASON_LIST_RED.append(season)
            SEASON_LIST = SEASON_LIST_RED
                 
        
        REQ_LIST=[]
        for season in SEASON_LIST:
            m = requestors.Season_req(season[0],season[1])
            REQ_LIST.append(m)
        return REQ_LIST
            
            
        SEASON_LIST=[]
        return SEASON_LIST
    
    def getOwnList(self):
        '''
        Not useful for time aggregation, but to get requestors in order to match with observations'''
        
        REQ_LIST=[]
        if self.inputFrequency == 'daily':
            for t in self.Timelist:
                REQ_LIST.append(requestors.Daily_req(t.year,t.month,t.day))
        if self.inputFrequency == 'monthly':
            for t in self.Timelist:
                REQ_LIST.append(requestors.Monthly_req(t.year, t.month))
        if self.inputFrequency == 'weekly' :
            for t in self.Timelist:
                REQ_LIST.append(requestors.Weekly_req(t.year, t.month, t.day)) 
        return REQ_LIST
    
    def getDecadalList(self):
        LIST=[]
        return LIST

    def couple_with(self,datetimelist):
        Coupled_List=[]
        for ir, req in enumerate(self.getOwnList()):
            LIST_of_IND=[]
            for ind_date, d in enumerate(datetimelist):
                if (d >= req.starttime) & (d <= req.endtime) :
                    LIST_of_IND.append(ind_date)
            if (len(LIST_of_IND) >0 ): Coupled_List.append((self.Timelist[ir],LIST_of_IND))
        return Coupled_List
        
                        
            
        



# sys.exit()
# TL = TimeList('20150101-00:00:00','20150701-00:00:00',dateformat, INPUTDIR,"ave*N1p.nc")
# M= TL.getWeekList('Tuesday')
# m=M[0]
# TL.select(m)
# 
# 
# for iSeas in range(4):
#     m = requestors.Clim_season('winter')
#     TL.select(m)
# 
# for imonth in range(12):
#     m = requestors.Clim_month(imonth)
#     TL.select(m)