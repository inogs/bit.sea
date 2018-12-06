import datetime

class TimeInterval():
    def __init__(self, starttime="19500101", endtime="21000101", dateformat='%Y%m%d'):
        self.start_time = datetime.datetime.strptime(starttime, dateformat)
        self.end_time = datetime.datetime.strptime(endtime   ,dateformat)
        assert self.start_time <= self.end_time
    
    def __repr__(self):
        return self.start_time.__repr__() + " , " + self.end_time.__repr__() 

    def contains(self, specific_time):
        '''
        Argument:
        * specific_time * : a datetime object
        Returns:
           True if specific_time is inside the time interval
        '''
        if (specific_time >= self.start_time) and (specific_time <= self.end_time):
            return True
        else:
            return False

    def overlapTime(self,T2):
        '''
          Argument:
          T2 : a TimeInterval object
          Returns:
            the number of seconds of the overlapping interval
        '''
    
    #------------------s1++++++++++e1--------------------- # Window of temporal average of model output file
    #                              | 
    #--------------------------s2++++++++++e2------------- # Window of aggregation defined in PLOT_LIST
    #                           |  |
    #--------------------------<++++>--------------------- # overlapping time window = theWindow
    
        theWindow = min(self.end_time,T2.end_time) - max(self.start_time,T2.start_time)     
    
    # Convert time window in seconds and nornamlize over the aggreation time window
        res       = float( theWindow.days*86400 + theWindow.seconds ) 
    
    # return only non-negative results
    
        return max( 0, res)
    
    def isInWindow(self,T2):
        '''
          Argument:
            T2 : a TimeInterval object
          Returns: True if there is an overlapping interval, False elsewhere.
        '''
        return self.overlapTime(T2) > 0
    
    def isInside(self,T2):
        raise NotImplementedError

    def length(self):
        '''
        Returns lengths of time interval in seconds
        '''
        theWindow = self.end_time - self.start_time
        res       = float( theWindow.days*86400 + theWindow.seconds ) 
        return res
    
    @staticmethod
    def fromdatetimes(timestart,time_end):
        '''
        Arguments:
        * timestart * a datetime object
        * time_end  * idem
        Returns a TimeInterval Object
        '''
        TI = TimeInterval()
        assert timestart <= time_end
        TI.start_time = timestart
        TI.end_time   = time_end
        return TI
    
    
