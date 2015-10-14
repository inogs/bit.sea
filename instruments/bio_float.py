import scipy.io.netcdf as NC
import numpy as np
import datetime
import os

from instruments.instrument import Instrument, Profile
from commons.time_interval import TimeInterval

class BioFloatProfile(Profile):
    def __init__(self, var, time, lon, lat, my_float, mean=None):
        self.var = var
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        self.mean = mean

        self._read_values = False
        self._pres = None
        self._data = None

    def _read_profile(self):
        self._pres, self._data = self._my_float.read(self.var, self.mean)
        self._read_values = True

    @property
    def pres(self):
        if not self._read_values:
            self._read_profile()
        return self._pres

    @property
    def data(self):
        if not self._read_values:
            self._read_profile()
        return self._data



class BioFloat(Instrument):

    default_mean = None

    def __init__(self,lon,lat,time,filename,available_params):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename
        self.available_params = available_params
        wmo, cycle = os.path.basename(filename).rsplit("_")
        self.wmo = wmo[2:]
        self.cycle = int(cycle[:3])

    def __searchVariable_on_parameters(self, var):
        '''
        returns index of profile on which variable has to be read,
        -1 if fails (PARAMETERS variable is blank)

        '''

        ncIN = NC.netcdf_file(self.filename,'r')
        N_PROF    = ncIN.dimensions['N_PROF']
        N_PARAM   = ncIN.dimensions['N_PARAM']
        PARAMETER = ncIN.variables['PARAMETER'].data.copy()
        ncIN.close()

        if PARAMETER.tostring().rstrip() == '':
            return -1
        else:
            for iprof in range(N_PROF):
                for iparam in range(N_PARAM):
                    s = PARAMETER[iprof,0,iparam,:].tostring().rstrip()
                    if s==var:
                        return iprof

    def __fillnan(self, ncObj,var):
        varObj = ncObj.variables[var]
        fillvalue = varObj._FillValue
        M = varObj.data.copy()
        M[M==fillvalue] = np.NaN;
        return M

    def __merge_profile_with_adjusted(self,profile, profile_adj):
        N_LEV = profile.size
        res = np.zeros_like(profile)
        if np.isnan(profile_adj).all():
            res = profile
        else:
            for k in range(N_LEV):
                if np.isnan(profile_adj[k]):
                    res[k] = profile[k]
                else:
                    res[k] = profile_adj[k]
        return res

    def __merge_var_with_adjusted(self, ncObj,var):
        N_PROF= ncObj.dimensions['N_PROF']
        M     = self.__fillnan(ncObj, var)
        M_ADJ = self.__fillnan(ncObj, var + "_ADJUSTED")
        M_RES = M
        for iprof in range(N_PROF):
            M_RES[iprof,:] = self.__merge_profile_with_adjusted(M[iprof,:], M_ADJ[iprof,:])

        return M_RES


    def read_raw(self,var):
        '''
        Reads data from file
        '''
        iProf = self.__searchVariable_on_parameters(var); #print iProf
        ncIN=NC.netcdf_file(self.filename,'r')

        if iProf== -1 :
            M = self.__merge_var_with_adjusted(ncIN, var)
            N_PROF= ncIN.dimensions['N_PROF']
            N_MEAS = np.zeros((N_PROF),np.int32)
            for iprof in range(N_PROF):
                N_MEAS[iprof] = (~np.isnan(M[iprof,:])).sum()
            #per il momento cerco il massimo
            iProf = N_MEAS.argmax()
            print "iprof new value", iProf

        M = self.__merge_var_with_adjusted(ncIN, var)
        PRES    = self.__merge_var_with_adjusted(ncIN, 'PRES')
        Profile =  M[iProf,:]
        Pres    =  PRES[iProf,:]
        ncIN.close()

        # Elimination of negative pressures or nans
        nanPres = np.isnan(Pres)
        Pres[nanPres] = 1 # just for not complaining
        badPres    = (Pres<=0) | nanPres
        badProfile = np.isnan(Profile)
        bad = badPres | badProfile

        return Pres[~bad], Profile[~bad]
        #return Pres, Profile

    def read(self, var, mean=None):
        '''
        Reads profile data from file and optionally applies a filter to the data
        '''
        pres, prof = self.read_raw(var)
        if mean == None:
            if BioFloat.default_mean != None:
                return pres, BioFloat.default_mean.compute(prof, pres)
            else:
                return pres, prof
        else:
            return pres, mean.compute(prof, pres)

    def plot(self,Pres,profile):
        pl.figure()
        pl.plot(profile,Pres)
        pl.gca().invert_yaxis()
        pl.show(block=False)

    def profiles(self, var, mean=None):
        return [BioFloatProfile(var, self.time, self.lon, self.lat, self, mean)]

    @staticmethod
    def from_file(filename):
        '''
        Returns the single Bio_Float instance corresponding to filename
        '''
        mydtype= np.dtype([
                  ('file_name','S200'),
                  ('lat',np.float32),
                  ('lon',np.float32),
                  ('time','S17'),
                  ('parameters','S200')] )

        FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Index.txt"
        INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
        nFiles=INDEX_FILE.size
        for iFile in range(nFiles):
            timestr          = INDEX_FILE['time'][iFile]
            lon              = INDEX_FILE['lon' ][iFile]
            lat              = INDEX_FILE['lat' ][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile]
            available_params = INDEX_FILE['parameters'][iFile]
            float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
            if filename == thefilename :
                        return Bio_Float(lon,lat,float_time,filename,available_params)
        return None


def FloatSelector(var, T, region):
    '''
    Arguments:
       var is a string indicating variable, 
          if var is None, no selection is done about variable
       T is as TimeInterval istance
       region is a region istance
    '''
    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )
    
    FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Index.txt"
    #FloatIndexer="Float_Index.txt"
    INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
    nFiles=INDEX_FILE.size    
    selected = []
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')

        if var is None :
            VarCondition = True
        else:
            VarCondition = var in available_params

        if VarCondition:
            if T.contains(float_time) and region.is_inside(lon, lat):
                selected.append(BioFloat(lon,lat,float_time,filename,available_params))
                      
    return selected
