import numpy as np
from instruments import var_conversions
class Profile(object):
    def read(self,var):
        raise NotImplementedError
    def name(self):
        raise NotImplementedError
    def ID(self):
        raise NotImplementedError
    def reference_var(self):
        raise NotImplementedError

class Instrument(object):
    pass



class ContainerProfile(Profile):
    def __init__(self,lon,lat,time, depth, values,name, datasetname):
        self.lon     = lon
        self.lat     = lat
        self.time    = time
        self.pres    = depth
        self.profile = values
        self._name    = name
        self.datasetname = datasetname
        self.has_adjusted = False

    def __eq__(self, other):
        if isinstance(other, ContainerProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                return True
        else:
            return False

    def read(self,var,read_adjusted=True):
        '''
        Here all arguments are unused. They are defined just to be compliant with other Profile objects
        (such biofloats ) where arguments are meaningful.

        Return pres, profile, Qc
        Qc is a dummy np.array of 2
        '''
        return self.pres, self.profile, np.ones_like(self.pres, np.int)*2

    def name(self):
        return self._name
    def ID(self):
        return  self._name + "_" + self.time.strftime("%Y%m%d_") + str(self.lon) + "_"+ str(self.lat)
    def reference_var(self, var):
        '''
        Returns the reference varname, for a given profile object and
        a ogstm model varname
        For BioFloats p.reference_var('O2o') returns 'DOXY'
        '''
        if (self.datasetname == 'Nutrients'): return var_conversions.NUTRVARS[var]
        if (self.datasetname == 'Carbon') : return var_conversions.CARBONVARS[var]
        if (self.datasetname == 'Massimili') : return var_conversions.MASSIMILIVARS[var]
        if (self.datasetname == 'Socat') : return var_conversions.SOCAT_VARS[var]

