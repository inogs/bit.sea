import datetime
from dateutil.relativedelta import relativedelta
from commons.time_interval import TimeInterval

class Clim_season():
    '''
    Requestor object  - for a climatologit season - used in Timelist.select() method.

    '''
    def __init__(self,season,seasonobj):
        '''
        Arguments:
        * season    * : integer
        * seasonobj * : an instance of class season, where features of season are defined
        
        Season object can be obtained as follows:
        import season
        seasonObj = season.season()
        
        Then seasonObj defines by default seasons like this:
        Winter = 0
        Spring = 1
        Summer = 2
        Fall   = 3
        
        Example:

        req = Clim_season(2,seasonObj)
        print(req.string)
        
        
        '''
        assert season in range(seasonobj.get_seasons_number())
        self.seasonobj = seasonobj
        self.season=season
        a=Season_req(2000,season,seasonobj)
        self.string = a.longname
    def __repr__(self):
        return "Climatologic Seasonal requestor object : "  + self.string

    def contains(self,time):
        '''
        Argument:
        * time * : a datetime object
        Returns:
           True if time is inside every time interval of the season
        '''
        yearly_date = datetime.datetime(self.seasonobj._reference_year, time.month, time.day, time.hour, time.minute, time.second)
        TI, nameseas = self.seasonobj.get_season_dates(self.season)
        return TI.contains(yearly_date)


class Clim_month():
    '''
    Requestor object - for a climatologic month - used in Timelist.select() method.

    Example:

    req = Clim_month(2)
    '''
    def __init__(self,month):
        assert month in range(1,13)
        self.month = month
        self.string = '%02d' %month
    def __repr__(self):
        return "Climatologic Monthly requestor object : "  + self.string
    def longname(self):
        a = datetime.datetime(2000,self.month,1)
        return a.strftime("Clim_%b")
    def contains(self,time):
        '''
        Argument:
        * time * : a datetime object
        Returns:
           True if time is inside every time interval of the month
        '''
        montlhy_req=Monthly_req(time.year, self.month)
        return montlhy_req.time_interval.contains(time)
         

class Clim_day():
    ''' Requestor object for - climatologic day -used in Timelist.select() method.
         '''
    def __init__(self,month,day):
        a=datetime.datetime(2001,month,day) #vede se giorno e mese sono corretti
        self.month = month
        self.day   = day
        self.string = a.strftime("%m%d")
    def __repr__(self):
        return "Climatologic Daily requestor object : "  + self.string

class Clim_Interval_req():
    ''' Useful for clima days averages
    Interval_req(1,25,days=10)
    '''
    def __init__(self,month,day, hour=12, days=10):
        self.month  = month
        self.day    = day
        self.hour   = hour
        centertime     = datetime.datetime(2001,self.month,self.day,self.hour)
        delta = datetime.timedelta(days)
        self.time_interval = TimeInterval.fromdatetimes(centertime-delta/2, centertime+delta/2)
        self.string  = centertime.strftime("Climatological Interval requestor %m%d")
        self.deltadays = days
    def __repr__(self):
        return "Interval requestor object: " + self.string + "  delta :  " + str(self.deltadays)  + " days"

class Monthly_req():
    ''' Requestor object - for specific month - used in Timelist.select() method.
        '''
    def __init__(self,year,month):
        self.year=year
        self.month=month
        self.string = "%d%02d" % (self.year, self.month)
        t = TimeInterval()
        t.start_time  = datetime.datetime(self.year, self.month,1)
        t.end_time    = t.start_time + relativedelta(months = 1)
        self.time_interval = t
    def __repr__(self):
        return "Monthly requestor object: " + self.string


class Yearly_req():
    '''Requestor object - for specific year - used in Timelist.select() method.
    '''
    def __init__(self,year):
        self.year=year
        t = TimeInterval()
        t.start_time = datetime.datetime(year  ,1,1)
        t.end_time   = datetime.datetime(year+1,1,1)
        self.time_interval =t
        self.string    = str(year)
    def __repr__(self):
        return "Yearly requestor object: " + self.string

class Weekly_req():
    ''' Requestor object - for a specific week - used in Timelist.select() method.
    The week is indicated by its central day.

    Example:

    req=Weekly_req(2015,3,5)
    print(r.isoweekday)
    '''
    def __init__(self,year,month,day):
        self.year   = year
        self.month  = month
        self.day    = day
        centertime     = datetime.datetime(self.year,self.month,self.day,12)
        t = TimeInterval()
        deltaseconds = 3.5*24*3600
        t.start_time = centertime - datetime.timedelta(seconds=deltaseconds)
        t.end_time   = centertime + datetime.timedelta(seconds=deltaseconds-1)
        self.time_interval = t
        self.string  = centertime.strftime("%Y%m%d")
        self.isoweekday = centertime.isoweekday()
    def __repr__(self):
        return "Weekly requestor object: " + self.string


class Daily_req():
    ''' Requestor object - for a specific day - used in Timelist.select() method.

    Example:

    req=Daily_req(2015,3,5)
    '''
    def __init__(self,year,month,day):
        self.year   = year
        self.month  = month
        self.day    = day
        t = TimeInterval()
        t.start_time =  datetime.datetime(self.year,self.month,self.day,0)
        t.end_time   = t.start_time + datetime.timedelta(days=1)
        self.time_interval = t
        self.string  = t.start_time.strftime("%Y%m%d")
    def __repr__(self):
        return "Daily requestor object: " + self.string

class Season_req():
    '''
    Requestor object - for a specific season - used in Timelist.select() method.

    '''
    def __init__(self,year,num_season,seasonobj):
        '''       
         Arguments:
        * year      * : integer
        * season    * : integer
        * seasonobj * : an instance of class season, where features of season are defined
        
        Season object can be obtained as follows:
        import season
        seasonObj = season.season()
        
        Then seasonObj defines by default seasons like this:
        Winter = 0
        Spring = 1
        Summer = 2
        Fall   = 3
        
        Example:

        req = Season_req(2012,2,seasonObj)
        print(req.longname)
        '''
        self.year   = year
        self.season = num_season
        ti_ref,self.longname = seasonobj.get_season_dates(num_season)

        delta_years=year - seasonobj._reference_year

        t = TimeInterval()
        t.start_time  = ti_ref.start_time + relativedelta(years=delta_years)
        t.end_time    = ti_ref.end_time   + relativedelta(years=delta_years)
        self.time_interval = t
        self.string   = str(year) + " " + self.longname[:3]


    def __repr__(self):
        return "Season requestor object :" + self.string






class Decadal_req():
    ''' Decadal requestor
    Decadal_req(2010) requests years 2010,2011,..., 2019
    Decadal_req(2011) requests years 2011,2012,..., 2020
    '''

    def __init__(self,decad):

        self.decad = decad

        # Decades are defined to start at 0 up to 9 (e.g.,2010-2019)
        # or start at 1 up to 10 (e.g., 2011-2020)
        # anything else is considered a generic interval.
        q,r=divmod(decad, 10)
        if not r in [0,1]:
            raise ValueError("Input is not a decade. Are you trying to implement an interval?")
        t = TimeInterval()
        self.startyear = decad
        self.end__year = self.startyear+9
        t.start_time = datetime.datetime(self.startyear  ,1,1,0,0,0)
        t.end_time   = datetime.datetime(self.end__year+1,1,1,0,0,0)
        self.time_interval = t
    def __repr__(self):
        return "Decadal requestor object : %d ... %d"  %(self.startyear, self.end__year)


class Interval_req():
    ''' Useful for 10 days averages
    Interval_req(2005,1,25,days=10)
    '''
    def __init__(self,year,month,day, hour=12, days=10):
        self.year   = year
        self.month  = month
        self.day    = day
        self.hour   = hour
        centertime     = datetime.datetime(self.year,self.month,self.day,self.hour)
        delta = datetime.timedelta(days)
        self.time_interval = TimeInterval.fromdatetimes(centertime-delta/2, centertime+delta/2)
        self.string  = centertime.strftime("%Y%m%d")
        self.deltadays = days
    def __repr__(self):
        return "Interval requestor object: " + self.string + "  delta :  " + str(self.deltadays)  + " days"


class Hourly_req():
    '''
    '''
    def __init__(self,year,month,day, hour, delta_hours=1):
        self.year   = year
        self.month  = month
        self.day    = day
        self.hour   = hour
        self.centertime     = datetime.datetime(self.year,self.month,self.day,self.hour)
        delta = datetime.timedelta(hours=delta_hours)

        self.time_interval = TimeInterval.fromdatetimes(self.centertime-delta/2, self.centertime+delta/2)
        self.string  = self.centertime.strftime("%Y%m%d-%H%M%S")
        self.delta_hours = delta_hours
    def __repr__(self):
        return "Hourly requestor object: " + self.string + "  delta :  " + str(self.delta_hours)  + " hours"

class seconds_req():
    '''
    '''
    def __init__(self,year,month,day, hour, minutes, delta_seconds=1800):
        self.year    = year
        self.month   = month
        self.day     = day
        self.hour    = hour
        self.minutes = minutes
        centertime    = datetime.datetime(self.year,self.month,self.day,self.hour,self.minutes)
        delta = datetime.timedelta(seconds=delta_seconds)

        self.time_interval = TimeInterval.fromdatetimes(centertime-delta/2, centertime+delta/2)
        self.string  = centertime.strftime("%Y%m%d-%H%M%S")
        self.delta_seconds = delta_seconds
    def __repr__(self):
        return "Seconds requestor object: " + self.string + "  delta :  " + str(self.delta_seconds)  + " seconds"



class Clim_Hourly_req():
    '''
    '''
    def __init__(self, hour, delta_hours=1):
        self.hour   = hour
        self.string  = ""
        self.delta_hours = delta_hours
    def __repr__(self):
        return "CLimatologic requestor object: " + self.string + "  delta :  " + str(self.delta_hours)  + " hours"





class Generic_req():
    ''' Generic requestor
        Based on timeinterval object, for non standard times.
        For example, an year defined between
    '''

    def __init__(self,ti):

        self.time_interval = ti
    def __repr__(self):
        return "Generic requestor object definded by: %s  "  %self.time_interval.__repr__()
    
if __name__=="__main__":
    a = Interval_req(2015,1,25,days=9)
