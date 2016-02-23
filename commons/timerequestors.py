import datetime
from dateutil.relativedelta import relativedelta
from time_interval import TimeInterval

class Clim_season():
    ''' 
    Requestor object  - for a climatologit season - used in Timelist.select() method.
    
    Assumption for integers indicating seasons:
    Winter = 0
    Spring = 1
    Summer = 2
    Fall   = 3
    
    Example:
    
    req= Clim_req(2) 
    print req.string
    ''' 
    def __init__(self,season):
        assert season in range(4)
        self.season=season
        a=Season_req(2000,season)
        self.string = a.longname
    def __repr__(self):
        return "Climatologic Seasonal requestor object : "  + self.string

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
    print r.isoweekday
    '''
    def __init__(self,year,month,day):
        self.year   = year
        self.month  = month
        self.day    = day
        centertime     = datetime.datetime(self.year,self.month,self.day,12)
        t = TimeInterval() 
        t.start_time = centertime - datetime.timedelta(days=3.5)
        t.end_time   = centertime + datetime.timedelta(days=3.5)
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
    
    Assumption for integers indicating seasons:
        Winter = 0
        Spring = 1
        Summer = 2
        Fall   = 3
    
    Example:
    
    req= Season_req(2012,2)
    print req.longname 
        ''' 
    def __init__(self,year,season):
        self.year   = year
        self.season = season
        
        assert season in [0,1,2,3]
        t = TimeInterval()
        if season == 0: #win
            t.start_time  = datetime.datetime(year ,  1, 1, 0, 0)
            t.end_time    = datetime.datetime(year ,  4, 1, 0, 0)
            self.string   = '%d win' %year
            self.longname = "Winter"
        if season == 1 : 
            t.start_time  = datetime.datetime(year ,  4, 1, 0, 0)
            t.end_time    = datetime.datetime(year ,  7, 1, 0, 0)
            self.string   = '%d spr' %year
            self.longname = 'Spring'
        if season == 2 :
            t.start_time  = datetime.datetime(year ,  7, 1, 0, 0)
            t.end_time    = datetime.datetime(year , 10, 1, 0, 0)
            self.string   = '%dsum' %year
            self.longname = 'Summer'
            
        if season == 3 :
            t.start_time  = datetime.datetime(year , 10, 1, 0, 0)
            t.end_time    = datetime.datetime(year+1, 1, 1, 0, 0)
            self.string   = '%d aut' %year
            self.longname = 'Autumn'
        
        self.timeinterval = t
            
    def __repr__(self):
        return "Season requestor object :" + self.string       
        

            
            
        

class Decadal_req():
    ''' Decadal requestor
    Decadal_req(2010) requests years 2010,2011,..., 2019
    Decadal_req(2011) requests years 2011,2012,..., 2020 
    '''
    
    def __init__(self,decad):

        self.decad = decad
        
        q,r=divmod(decad, 10)
        assert r in [0,1]
        t = TimeInterval()
        self.startyear = decad                   
        self.end__year = self.startyear+9
        t.start_time = datetime.datetime(self.startyear  ,1,1,0,0,0)
        t.end_time   = datetime.datetime(self.end__year+1,1,1,0,0,0)
        self.timeinterval = t
    def __repr__(self):
        return "Decadal requestor object : %d ... %d"  %(self.startyear, self.end__year)


class Interval_req():
    ''' Useful for 10 days averages
    Interval_req(2005,1,25,'days=10')
    '''
    def __init__(self,year,month,day, deltastr):
        self.year   = year
        self.month  = month
        self.day    = day
        centertime     = datetime.datetime(self.year,self.month,self.day,12)
        delta = relativedelta(10)
        exec 'delta= relativedelta(' + deltastr + ')'
        self.timeinterval = TimeInterval.fromdatetimes(centertime-delta/2, centertime+delta/2)
        self.string  = centertime.strftime("%Y%m%d")
        self.deltastr = deltastr
    def __repr__(self):
        return "Interval requestor object: " + self.string + "  delta :  " + self.deltastr 




class Generic_req():
    ''' Generic requestor
        Based on timeinterval object, for non standard times.
        For example, an year defined between 
    '''
    
    def __init__(self,ti):

        self.timeinterval = ti
    def __repr__(self):
        return "Generic requestor object definded by: %s  "  %self.timeinterval.__repr__() 




    