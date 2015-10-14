from instruments.instrument import Instrument, Profile

class MooringProfile(Profile):
    def __init__(self,time,lat,lon,filename):
        self.time = time
        self.lat  = lat
        self.lon  = lon
        self.filename = filename

class Mooring(Instrument):
    def __init__(self, lon, lat, name):
        self.lon = lat
        self.lat = lat
        self.name = name
    
    def profiles(self):
        raise NotImplementedError
    
    def read_raw(self, var, profile=None):
        '''
        If profile is None the entire 2d array (depth,time) is returned
        If a profile object is provided, a 1d array (dimension depth) is returned
        corresponding to the pointprofile provided
        '''
        raise NotImplementedError
        if profile is None:
            pass
            
        #returns array(times,depths) 
    
    def read(self, var, mean=None, profile=None):
        raise NotImplementedError
        #returns array(times,depths) 


