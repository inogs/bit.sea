class Profile(object):
    def read(self,var):
        raise NotImplementedError
    def name(self):
        raise NotImplementedError

class Instrument(object):
    pass



class ContainerProfile(Profile):
    def __init__(self,lon,lat,time, depth, values,name):
        self.lon     = lon
        self.lat     = lat
        self.time    = time
        self.pres    = depth
        self.profile = values
        self.name    = name

    def read(self,var):
        return self.pres, self.profile

