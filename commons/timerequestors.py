import datetime
from dateutil.relativedelta import relativedelta


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
    def __str__(self):
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
    def __str__(self):
        return "Climatologic Monthly requestor object : "  + self.string  

class Clim_day():
    ''' Requestor object for - climatologic day -used in Timelist.select() method.
         '''    
    def __init__(self,month,day):
        a=datetime.datetime(2001,month,day) #vede se giorno e mese sono corretti
        self.month = month
        self.day   = day
        self.string = a.strftime("%m%d")
    def __str__(self):
        return "Climatologic Daily requestor object : "  + self.string        

class Monthly_req():
    ''' Requestor object - for specific month - used in Timelist.select() method.
        '''
    def __init__(self,year,month):
        self.year=year
        self.month=month
        self.string = "%d%02d" % (self.year, self.month)
        self.starttime  = datetime.datetime(self.year, self.month,1)
        self.endtime    = self.starttime + relativedelta(months = 1)
    def __str__(self):
        return "Monthly requestor object: " + self.string   
        
        
class Yearly_req():
    '''Requestor object - for specific year - used in Timelist.select() method.
    '''
    def __init__(self,year):
        self.year=year
        self.starttime = datetime.datetime(year  ,1,1)
        self.endtime   = datetime.datetime(year+1,1,1)
        self.string    = str(year)
    def __str__(self):
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
        self.starttime = centertime - datetime.timedelta(days=3.5)
        self.endtime = centertime + datetime.timedelta(days=3.5)
        self.string  = centertime.strftime("%Y%m%d")
        self.isoweekday = centertime.isoweekday()
    def __str__(self):
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
        self.starttime =  datetime.datetime(self.year,self.month,self.day,0)
        self.endtime = self.starttime + datetime.timedelta(days=1)
        self.string  = self.starttime.strftime("%Y%m%d")
    def __str__(self):
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
        
        if season == 0: #win
            self.starttime = datetime.datetime(year ,  1, 1, 0, 0)
            self.endtime   = datetime.datetime(year ,  4, 1, 0, 0)
            self.string    = 'win %d' %year 
            self.longname  = "Winter"
        if season == 1 : 
            self.starttime = datetime.datetime(year ,  4, 1, 0, 0)
            self.endtime   = datetime.datetime(year ,  7, 1, 0, 0)
            self.string    = 'spr %d' %year
            self.longname  = 'Spring'
        if season == 2 :
            self.starttime = datetime.datetime(year ,  7, 1, 0, 0)
            self.endtime   = datetime.datetime(year , 10, 1, 0, 0)
            self.string    = 'sum %d' %year
            self.longname  = 'Summer'
            
        if season == 3 :
            self.starttime = datetime.datetime(year , 10, 1, 0, 0)
            self.endtime   = datetime.datetime(year+1, 1, 1, 0, 0)
            self.string    = 'fal %d' %year 
            self.longname  = 'Fall'
            
    def __str__(self):
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
        
        self.startyear = decad                   
        self.end__year = self.startyear+9
        self.starttime = datetime.datetime(self.startyear  ,1,1,0,0,0)
        self.endtime   = datetime.datetime(self.end__year+1,1,1,0,0,0)
    def __str__(self):
        return "Decadal requestor object : %d ... %d"  %(self.startyear, self.end__year)






    