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
        self._name    = name

    def __eq__(self, other):
        if isinstance(other, ContainerProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                True
        else:
            return False

    def read(self,var):
        return self.pres, self.profile

    def name(self):
        return self._name

