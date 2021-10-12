from __future__ import print_function
from commons import timerequestors as requestors
from commons import genUserDateList as DL
import os,glob
import datetime
import numpy as np
from commons import season
from commons import IOnames
from commons.time_interval import TimeInterval
from commons.utils import addsep

seasonobj = season.season()

def computeTimeWindow(freqString,currentDate):

    if (freqString == 'daily'):   req = requestors.Daily_req(currentDate.year,currentDate.month,currentDate.day)
    if (freqString == 'weekly'):  req = requestors.Weekly_req(currentDate.year, currentDate.month,currentDate.day)
    if (freqString == 'monthly'): req = requestors.Monthly_req(currentDate.year, currentDate.month)
    if (freqString == 'yearly'):  req = requestors.Yearly_req(currentDate.year)
    if freqString.startswith('days'):
        pos=freqString.find("=")
        ndays=int(freqString[pos+1:])
        req = requestors.Interval_req(currentDate.year,currentDate.month,currentDate.day,days=ndays)
    if freqString.startswith('hours'):
        pos=freqString.find("=")
        nhours=int(freqString[pos+1:])
        req = requestors.Hourly_req(currentDate.year,currentDate.month,currentDate.day,currentDate.hour, delta_hours=nhours)
    if freqString.startswith('seconds='):
        pos=freqString.find("=")
        nseconds=int(freqString[pos+1:])
        req = requestors.seconds_req(currentDate.year,currentDate.month,currentDate.day,currentDate.hour,currentDate.minute, delta_seconds=nseconds) 
    return TimeInterval.fromdatetimes(req.time_interval.start_time, req.time_interval.end_time)

class TimeList():

    def __init__(self,datelist,forceFrequency=None):
        '''
        TimeList object is created by providing a list of datetime objects
        (At least 2).
        '''
        nTimes = len(datelist)
        self.Timelist = datelist
        self.Timelist.sort()
        self.nTimes   = nTimes

        self.inputdir     = None
        self.searchstring = None
        self.filelist     = None
        self.filtervar    = None
        self.inputFrequency = None
        if forceFrequency is not None:
            self.inputFrequency = forceFrequency
        else:
            if (nTimes > 1 ) :
                self.inputFrequency= self.__searchFrequency()
                self.timeinterval = TimeInterval.fromdatetimes(self.Timelist[0], self.Timelist[-1])

    @staticmethod
    def fromfilenames(timeinterval, inputdir,searchstring, filtervar=None, prefix='ave.', dateformat="%Y%m%d-%H:%M:%S",hour=12,forceFrequency=None):
        '''
        Generates a TimeList object by reading a directory

        HYPOTHESIS:
         - the Input directory has files containing a date in their names
         - every file refers to a single time frame
         - The date of in the file name can be considered as centered in their period
        The generated datetime list has all the files concerning the period indicated in the timeinterval.
        Some of these files can have the centered date out of that period.

        Arguments:
        * timeinterval   * a TimeInterval object
                           if it's None, no selection of time will be perfomed, 
                           all files in that directory will be taken in account.
        * inputdir       * string, directory to where files are
        * searchstring   * string, used to filter vars
        * filtervar      * string, e.g. 'N1p', used for filter a filelist
                           a file is vaild if the filename string contains filtervar
        * prefix         * string, the part of the filename before the date
        * dateformat     * string
        * hour           * integer
        * forceFrequency * string, like 'daily','weekly','monthly'. The default is None.
                           If set to None, the frequency is automatically calculated, else is forced to
                           the string provided. This is useful when the timeseries has gaps and is not continue and automatic calculation
                           can fail.


        Example:

        INPUTDIR="/pico/scratch/userexternal/gbolzon0/Carbonatic/wrkdir/MODEL/AVE_FREQ_1/"
        Time_int = TimeInterval('20141201-00:00:00','20150701-00:00:00',"%Y%m%d-%H:%M:%S")
        TL = TimeList.fromfilenames(Time_int, INPUTDIR,"ave*N1p.nc")
        TL = TimeList.fromfilenames(Time_int, INPUTDIR,"ave*nc",filtervar="N1p")

        For Sat data
        TL = TimeList.fromfilenames(None, INPUTDIR,"*nc", prefix='',dateformat='%Y%m%d')

        For physical forcings
        TL = TimeList.fromfilenames(None  , INPUTDIR,"T*.nc",prefix='T',dateformat='%Y%m%d')
        '''
        if timeinterval is None:
            timeinterval=TimeInterval("1900","2200","%Y")
        IOname = IOnames.filenamer(prefix,dateformat)
        if not os.path.exists(inputdir):
            raise NameError("Not existing directory " + inputdir)

        inputdir = addsep(inputdir)
        filelist_ALL = glob.glob(inputdir + searchstring)
        if not filtervar is None:
            filename, file_extension = os.path.splitext(filelist_ALL[0])
            #filelist_ALL=[f for f in filelist_ALL if f.endswith("." + filtervar + file_extension) ]
            filelist_ALL=[f for f in filelist_ALL if filtervar+".nc" in os.path.basename(f) ]
        assert len(filelist_ALL) > 0
        filenamelist=[]
        datetimelist=[]
        External_filelist=[]
        External_timelist=[]
        for pathfile in filelist_ALL:
            filename   = os.path.basename(pathfile)
            datestr     = filename[IOname.date_startpos:IOname.date_endpos]
            try:
                actualtime = datetime.datetime.strptime(datestr,IOname.dateformat)
                if timeinterval.contains(actualtime) :
                    filenamelist.append(pathfile)
                    datetimelist.append(actualtime)
                else:
                    External_filelist.append(pathfile)
                    External_timelist.append(actualtime)
            except:
                print("Warning: " + datestr + " does not exist!")

        TimeListObj = TimeList(datetimelist,forceFrequency=forceFrequency)
        filenamelist.sort()
        TimeListObj.timeinterval = timeinterval
        TimeListObj.inputdir     = inputdir
        TimeListObj.searchstring = searchstring
        TimeListObj.prefix       = prefix
        TimeListObj.filtervar    = filtervar
        TimeListObj.filelist = filenamelist



        if TimeListObj.inputFrequency is not None:
            for iFrame, t in enumerate(External_timelist):
                TimInt = computeTimeWindow(TimeListObj.inputFrequency, t)
                if TimeListObj.timeinterval.isInWindow(TimInt):
                    TimeListObj.filelist.append(External_filelist[iFrame])
                    TimeListObj.Timelist.append(External_timelist[iFrame])
            TimeListObj.filelist.sort()
            TimeListObj.Timelist.sort()
        TimeListObj.nTimes   =len(TimeListObj.filelist)


        # we force daily datetimes to have hours = 12
        if TimeListObj.inputFrequency == 'daily':
            for iFrame, t in enumerate(TimeListObj.Timelist):
                if t.hour==0:
                    newt = datetime.datetime(t.year,t.month,t.day,hour,0,0)
                    TimeListObj.Timelist[iFrame] = newt

        return TimeListObj


    def __searchFrequency(self):
        '''
        Returns strings: 'daily','weekly','monthly','hourly','10days'
        '''
        if len(self.Timelist)<2:
            timestr = self.timeinterval.start_time.strftime(" between %Y%m%d and ") +  self.timeinterval.end_time.strftime("%Y%m%d")
            print("Frequency cannot be calculated in " + self.inputdir + timestr)
            return None

        DIFFS = np.zeros((self.nTimes-1),np.float64)
        for iTime in range(1,self.nTimes):
            mydiff = self.Timelist[iTime]-self.Timelist[iTime-1]
            DIFFS[iTime-1] = mydiff.total_seconds()

        UniqueDiffs, nRepetions = np.unique(DIFFS,return_counts=True)
        moda = UniqueDiffs[nRepetions.argmax()]
        days = moda/(3600.*24)
        if days == 1:
            return "daily"
        if (days > 6 ) & (days < 8 ) :
            return "weekly" #centered in " + str(self.Timelist[1].isoweekday())
        if (days > 26) & (days < 32) :
            return "monthly"
        if (days < 1 ) & ( days >= 1./24.):
            return "hours=%d" %int(days*24)
        if (days < 1./24 ) :
            return "seconds=%d" % int(days*24.*3600.)
        if (days > 364 ) & (days < 367): 
            return "yearly"
        if (abs (int(days) -days)) < 0.1:
            return "days=%d" %days
        if days == 10:
            return "10days"
        if (days>1) & (days<7):
            return None
        else:
            raise NotImplementedError
            #hours = mydiff.seconds/3600
            #we want an integer number of hours
            #if (float(mydiff.seconds)/3600. == hours):


    def __generaltimeselector(self,requestor):
            SELECTION=[]
            weights  =[]

            if self.inputFrequency == "daily":
                for it, t in enumerate(self.Timelist):
                    if requestor.time_interval.contains(t):
                        SELECTION.append(it)
                        weights.append(1.)
            if self.inputFrequency in ['weekly','monthly','yearly','10days','days=3']:
                for it, t in enumerate(self.Timelist):
                    t1 = computeTimeWindow(self.inputFrequency,t)
                    t2 = TimeInterval.fromdatetimes(requestor.time_interval.start_time, requestor.time_interval.end_time)
                    weight = t1.overlapTime(t2)
                    if (weight > 0. ) :
                        SELECTION.append(it)
                        weights.append(weight)
            return SELECTION , np.array(weights)

    def select_one(self,requestor):
        '''
        Used to select a single time (or file) regardless to time aggregation

        index = select_one(requestor)
        Returned values:
          - an integer index indicating to access selected times (or files)


        '''
        assert not isinstance(requestor, requestors.Clim_day)
        assert not isinstance(requestor, requestors.Clim_month)
        assert not isinstance(requestor, requestors.Clim_season)
        for it, t in enumerate(self.Timelist):
            if requestor.time_interval.contains(t):
                return it
        print("Time not found")
        return None

    def selectWeeklyGaussWeights(self,requestor,std):
        '''
        Method for time aggregation for weekly req 
        with weigths from Gaussian distribution centered on week central date
        shape of Gaussian distribution from std input
        indexes, weights = select(requestor,std)
        Returned values:
         - a list of indexes (integers) indicating to access selected times (or files)
         - a numpy array of weights


        '''


        if isinstance(requestor, requestors.Weekly_req):
            assert self.inputFrequency != "monthly"
            assert self.inputFrequency != "weekly"
            assert self.inputFrequency != "10days"

            from scipy.stats import norm

            SELECTION=[]
            weights  =[]

            gaussweight = [norm.pdf(x,0,std) for x in range(4)]

            if self.inputFrequency in ["daily","days=2"]:
                for it,t in enumerate(self.Timelist):
                    if requestor.time_interval.contains(t):
                        SELECTION.append(it)
                        dayd = abs((t-requestor.time_interval.start_time).days-3)
                        weights.append(gaussweight[dayd])
            return SELECTION , np.array(weights)


        else: 
                raise NotImplementedError

    def select(self,requestor):
        '''
        Method for time aggregation
        indexes, weights = select(requestor)
        Returned values:
         - a list of indexes (integers) indicating to access selected times (or files)
         - a numpy array of weights


        '''
        if isinstance(requestor,requestors.seconds_req):
            assert self.inputFrequency in ["seconds=900", "seconds=1800"]
            SELECTION=[]
            weights = []
            for it,t in enumerate(self.Timelist):
                if requestor.time_interval.contains(t):
                    SELECTION.append(it)
                    weights.append(1.)
            return SELECTION , np.array(weights)

        if isinstance(requestor,requestors.Hourly_req):
            assert self.inputFrequency in [ 'hourly', "seconds=900", "seconds=1800"]
            SELECTION=[]
            weights = []
            for it,t in enumerate(self.Timelist):
                if requestor.time_interval.contains(t):
                    SELECTION.append(it)
                    weights.append(1.)
            return SELECTION , np.array(weights)

        if isinstance(requestor,requestors.Daily_req):
            # hourly values are treated as instantaneous values, not time averages
            assert self.inputFrequency in ["hourly","daily", "seconds=900", "seconds=1800"] # it does not matter how many hours
            SELECTION=[]
            weights = []
            for it,t in enumerate(self.Timelist):
                if requestor.time_interval.contains(t):
                    SELECTION.append(it)
                    weights.append(1.)
            return SELECTION , np.array(weights)


        if isinstance(requestor, requestors.Monthly_req):
            SELECTION=[]
            weights = []

            if self.inputFrequency in  ['daily','hourly',"seconds=900", "seconds=1800"] :
                for it, t in enumerate(self.Timelist):
                    if (t.year==requestor.year) & (t.month==requestor.month):
                        SELECTION.append(it)
                        weights.append(1.)
                return SELECTION , np.array(weights)
            if self.inputFrequency == 'weekly':
                for it,t in enumerate(self.Timelist):
                    t1 = computeTimeWindow("weekly",t);
                    t2 = TimeInterval.fromdatetimes(requestor.time_interval.start_time, requestor.time_interval.end_time)
                    weight = t1.overlapTime(t2);
                    if (weight > 0. ) :
                        SELECTION.append(it)
                        weights.append(weight)
                return SELECTION , np.array(weights)
            if self.inputFrequency[:5]=='days=':
                for it,t in enumerate(self.Timelist):
                    t1 = computeTimeWindow(self.inputFrequency,t);
                    t2 = TimeInterval.fromdatetimes(requestor.time_interval.start_time, requestor.time_interval.end_time)
                    weight = t1.overlapTime(t2);
                    if (weight > 0. ) :
                        SELECTION.append(it)
                        weights.append(weight)
                return SELECTION , np.array(weights)
            if self.inputFrequency == 'monthly':
                for it,t in enumerate(self.Timelist):
                    if requestor.time_interval.contains(t):
                        SELECTION.append(it)
                        weights.append(1.)
                return SELECTION , np.array(weights)
            else:
                raise NotImplementedError


        if isinstance(requestor, requestors.Weekly_req):
            assert self.inputFrequency != "monthly"
            assert self.inputFrequency != "weekly"
            assert self.inputFrequency != "10days"

            SELECTION=[]
            weights  =[]

            if self.inputFrequency in ["daily","days=2"]:
                for it,t in enumerate(self.Timelist):
                    if requestor.time_interval.contains(t):
                        SELECTION.append(it)
                        weights.append(1.)
            return SELECTION , np.array(weights)

        if isinstance(requestor, requestors.Interval_req):
            return self.__generaltimeselector(requestor)

        if isinstance(requestor, requestors.Season_req):
            return self.__generaltimeselector(requestor)

        if isinstance(requestor, requestors.Yearly_req):
            return self.__generaltimeselector(requestor)

        if isinstance(requestor, requestors.Decadal_req):
            return self.__generaltimeselector(requestor)
        if isinstance(requestor, requestors.Generic_req):
            return self.__generaltimeselector(requestor)

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

        if isinstance(requestor,requestors.Clim_Interval_req):
            SELECTION = []
            weights   = []
            YEARLIST=self.getYearlist()
            for year_req in YEARLIST:
                req = requestors.Interval_req(year_req.year, \
                                              requestor.month, \
                                              requestor.day, \
                                              requestor.deltadays)
                s,w = self.__generaltimeselector(req)
                SELECTION.extend(s)
                weights.extend(w)

            return SELECTION , np.array(weights)

        if isinstance(requestor,requestors.Clim_season):
            SELECTION = []
            weights   = []
            YEARLIST=self.getYearlist()
            for year_req in YEARLIST:
                req = requestors.Season_req(year_req.year, requestor.season,requestor.seasonobj)
                s,w = self.select(req)
                SELECTION.extend(s)
                weights.extend(w)
            return SELECTION , np.array(weights)
        if isinstance(requestor, requestors.Clim_Hourly_req):
            SELECTION=[]
            weights  =[]
            DAILYLIST=self.getDailyList()
            for day_req in DAILYLIST:
                req = requestors.Hourly_req(day_req.year, day_req.month, day_req.day, requestor.hour, requestor.delta_hours)
                s,w = self.select(req)
                SELECTION.extend(s)
                weights.extend(w)
            return SELECTION , np.array(weights)

    def getDailyList(self):
        '''
        Tested only for mooring case, interval = 3 hours
        '''
        REQ_LIST=[]
        t = self.Timelist[0]
        starting_centered_day = datetime.datetime(t.year,t.month,t.day)
        TL=DL.getTimeList(starting_centered_day,self.Timelist[-1] , days=1)
        for t in TL:
            d = requestors.Daily_req(t.year,t.month,t.day)
            indexes,_ = self.select(d)
            if len(indexes)>0:
                REQ_LIST.append(d)
        return REQ_LIST

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
            if interested_weekday == self.timeinterval.start_time.isoweekday():
                index = counter

        starting_centered_day = self.timeinterval.start_time + datetime.timedelta(days=PossibleShifts[index])

        TL=DL.getTimeList(starting_centered_day,self.Timelist[-1] , days=7)
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
            firstMonth = datetime.datetime(self.timeinterval.start_time.year, self.timeinterval.start_time.month,1)
            lastMonth  = datetime.datetime(self.timeinterval.end_time.year,  self.timeinterval.end_time.month  ,1)
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
        YEARLIST=[self.timeinterval.start_time.year]
        for t in self.Timelist:
            if t.year not in YEARLIST: YEARLIST.append(t.year)
        REQ_LIST=[]
        for year in YEARLIST:
            m = requestors.Yearly_req(year)
            REQ_LIST.append(m)
        return REQ_LIST



    def getSeasonList(self,seasonobj=seasonobj, extrap=False):
        '''Returns an ordered list of requestors, Season_req objects.
       By setting extrap=True, this method extrapolates out of the indicated period (Starttime,EndTime)
       Example: if the input is weekly, centered in 20120301, and Starttime=20120301,
       then extrapolation will return also 201202, because of the part of the week before the centered time.
        '''
        t=self.Timelist[0]
        SEASON_LIST=[(t.year, seasonobj.findseason(t))]
        for t in self.Timelist:
            newSeason = (t.year, seasonobj.findseason(t))
            if newSeason not in SEASON_LIST: SEASON_LIST.append(newSeason)
        if extrap:
            pass
        else:
            SEASON_LIST_RED=[]
            firstSeason=requestors.Season_req(self.timeinterval.start_time.year, seasonobj.findseason(self.timeinterval.start_time),seasonobj)
            lastSeason = requestors.Season_req(self.timeinterval.end_time.year,  seasonobj.findseason(self.timeinterval.end_time),seasonobj)
            for season in SEASON_LIST:
                req = requestors.Season_req(season[0],season[1],seasonobj)

                if (req.time_interval.start_time >= firstSeason.time_interval.start_time) & (req.time_interval.end_time <=lastSeason.time_interval.end_time):
                    SEASON_LIST_RED.append(season)
            SEASON_LIST = SEASON_LIST_RED


        REQ_LIST=[]
        for season in SEASON_LIST:
            m = requestors.Season_req(season[0],season[1],seasonobj)
            REQ_LIST.append(m)
        return REQ_LIST


        SEASON_LIST=[]
        return SEASON_LIST

    def getOwnList(self):
        '''
        Not useful for time aggregation, but to get requestors in order to match with observations'''

        REQ_LIST=[]
        if self.inputFrequency.startswith('seconds='):
            pos=self.inputFrequency.find("=")
            nSeconds=int(self.inputFrequency[pos+1:])
            for t in self.Timelist:
                REQ_LIST.append(requestors.seconds_req(t.year,t.month,t.day,t.hour,t.minute,nSeconds))
            return REQ_LIST
        if self.inputFrequency.startswith('hours'):
            pos=self.inputFrequency.find("=")
            nHours=int(self.inputFrequency[pos+1:])
            for t in self.Timelist:
                REQ_LIST.append(requestors.Hourly_req(t.year,t.month,t.day,t.hour, delta_hours=nHours))
            return REQ_LIST
        if self.inputFrequency == 'daily':
            for t in self.Timelist:
                REQ_LIST.append(requestors.Daily_req(t.year,t.month,t.day))
            return REQ_LIST
        if self.inputFrequency == '10days':
            for t in self.Timelist:
                REQ_LIST.append(requestors.Interval_req(t.year, t.month, t.day, 'days=10'))
            return REQ_LIST
        if self.inputFrequency == 'monthly':
            for t in self.Timelist:
                REQ_LIST.append(requestors.Monthly_req(t.year, t.month))
            return REQ_LIST
        if self.inputFrequency == 'weekly' :
            for t in self.Timelist:
                REQ_LIST.append(requestors.Weekly_req(t.year, t.month, t.day))
            return REQ_LIST
        if self.inputFrequency == 'yearly' :
            for t in self.Timelist:
                REQ_LIST.append(requestors.Yearly_req(t.year))
            return REQ_LIST
#       if self.inputFrequency == None:
        raise NotImplementedError
#       return REQ_LIST


    def getSpecificIntervalList(self,deltadays=10,starttime="19971001-12:00:00"):
        '''
        Useful in case of 10 days average, for example
        '''
        deltastr = 'days=' + str(int(deltadays)) 
        REQ_LIST=[]
        dl=DL.getTimeList(starttime, self.timeinterval.end_time.strftime("%Y%m%d-%H:%M:%S"), deltastr)
        for dateobj in dl:
            req= requestors.Interval_req(dateobj.year,dateobj.month,dateobj.day, deltadays)
            REQ_LIST.append(req)
        return REQ_LIST


    def getDecadalList(self):
        LIST=[]
        raise NotImplementedError
        return LIST

    def couple_with(self,datetimelist):
        '''
        Performs the association between two timelists:

        - a regular one, usual result of model
        - an irregular one, as provided from measurements

        Input:
        * datetimelist * a list of datetime objects, usually the irregular times of profiles objects

        Returns:
        * Coupled_List * a list of tuples (datetime object, list_of_indices)
        where
        - the datetime object is an element of self.Timelist (model regular times) concerned by datetimelist
        - list_of_indices is a list of integers such that :
            datetimelist[list_of_indeces] belong to the datetime object corresponding timeinterval
        '''
        Coupled_List=[]
        for ir, req in enumerate(self.getOwnList()):
            LIST_of_IND=[]
            for ind_date, d in enumerate(datetimelist):
                if req.time_interval.contains(d):
                    LIST_of_IND.append(ind_date)
            if (len(LIST_of_IND) >0 ): Coupled_List.append((self.Timelist[ir],LIST_of_IND))
        return Coupled_List

    def merge_with(self,datetimelist):
        '''
        Merge with a datetimelist
        Returns a datetimelist with all the dates
         from both lists
         without repetitions
        '''
        Merged_List=[]
        Merged_List.extend(self.Timelist)
        Merged_List.extend(datetimelist)
        Merged_List.sort()
        Merged_List=list(np.unique(Merged_List))
        return Merged_List

    def find(self,datetimeObj,returndiff=False):
        '''
        Finds the nearest
        Argument:
         a datetime object
        Returns the index of the nearest element
        '''
        D=np.zeros(self.nTimes)
        for i, d in enumerate(self.Timelist):
            diff=d-datetimeObj
            D[i]=np.abs(diff.total_seconds())
        if returndiff:
            return D.argmin(),D[D.argmin()]
        else:
            return D.argmin()



if __name__ == '__main__':
    Days17 = DL.getTimeList("19970127-00:00:00","19970601-00:00:00", days=17)
    TTL     = TimeList(Days17)
    MyReqList = TTL.getMonthlist()
    for req in MyReqList:
        ii,weights = TTL.select(req)

    H2 = DL.getTimeList("20180301-00:00:00","20200310-00:00:00", hours=2)
    TL     = TimeList(H2)
    REQS = TL.getOwnList()
    Min15 = DL.getTimeList("20180301-00:00:00","20200310-00:00:00", minutes=15)
    TL     = TimeList(Min15)
    Sec_req=requestors.seconds_req(2018,3,5,12,0,delta_seconds=900)
    ii,w = TL.select(Sec_req)

    Min15 = DL.getTimeList("20180301-00:00:00","20200310-00:00:00", minutes=15)
    TL     = TimeList(Min15)
    Hourly_req=requestors.Hourly_req(2018,3,5,12,delta_hours=2)
    ii,w = TL.select(Hourly_req)


    Clim_Hour_req  = requestors.Clim_Hourly_req(12,delta_hours=2)
    Clim_Montlhy_req = requestors.Clim_month(4)


    ii,w = TL.select(Clim_Hour_req)
    jj,wj = TL.select(Clim_Montlhy_req)
    intersection=[ k for k in ii if k in jj]


    TenDays = DL.getTimeList("19970127-00:00:00","19970601-00:00:00", days=10)
    TTL     = TimeList(TenDays)
    MyReqList = TTL.getMonthlist()
    for req in MyReqList:
        ii,weights = TTL.select(req)
        # print("\nii = ", ii, "\nw = ", weights)

    yearly=DL.getTimeList("19970101-00:00:00", "20150502-12:00:00", years=1)
    TLY = TimeList(yearly)
    REQSY=TLY.getOwnList()
    r=REQSY[0]
    ii,weights = TLY.select(r)

    monthly=DL.getTimeList("19970601-00:00:00", "20150502-12:00:00", months=1)
    TLM = TimeList(monthly)
    for iSeas in range(4):
        m = requestors.Clim_season(iSeas,seasonobj)
        TLM.select(m)

    for imonth in range(1,13):
        m = requestors.Clim_month(imonth)
        TLM.select(m)

    daily=DL.getTimeList("19970601-00:00:00", "20150502-12:00:00", days=1)
    TL = TimeList(daily)
    REQS=TL.getOwnList()
    r=REQS[0]
    ii,weights = TL.select(r)


    M= TL.getWeeklyList(2)
    m=M[0]
    TL.select(m)
    TL.getSeasonList(seasonobj)

    for iSeas in range(4):
        m = requestors.Clim_season(iSeas,seasonobj)
        TL.select(m)

    for imonth in range(1,13):
        m = requestors.Clim_month(imonth)
        TL.select(m)
    req = requestors.Monthly_req(2012,1)
    TLY.select_one(req)

