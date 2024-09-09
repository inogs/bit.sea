import numpy as np
import os,datetime
from instrument import Instrument, Profile
import scipy.io.netcdf as NC
import index_reader
import seawater
from seawater.library import T90conv


from commons.time_interval import TimeInterval

basetime = datetime.datetime(1950,1,1,0,0,0)
INDEX_FILE=index_reader.index_reader()
fillValue = -9999
#DOX1:units = "ml/l" ;
#CPHL:units = "milligram/m3" ;

class MooringProfile(Profile):
    def __init__(self,time,lon,lat, my_mooring,available_params):
        self.time = time
        self.lat  = lat
        self.lon  = lon
        self.available_params = available_params
        self._my_mooring = my_mooring

    def __eq__(self, other):
        if isinstance(other, MooringProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                return self._my_mooring == other._my_mooring
            else:
                return False
        else:
            return False

    def read(self,var, read_adjusted):
        '''
        read_adjusted is an unused parameter
        '''
        return self._my_mooring.read(var,self.time)

    def name(self):
        '''returns a string, the MyOcean identifier of the mooring.
        '''
        return self._my_mooring.name

class Mooring(Instrument):
    def __init__(self, lon, lat, filename, available_params):
        self.lon = lon
        self.lat = lat
        self.filename = filename
        self.available_params = available_params
        self.name = os.path.basename(filename)[10:-3]
        self.__file_already_read = False
        self._VAR = None
        self._QC  = None
        self._PRES = None
        self._timesInFile = None
        self._TEMP = None
        self._PSAL = None
        self._RHO  = None
        self._good_data = None

    def __eq__(self,other):
        if isinstance(other, Mooring):
            if (self.filename  == other.filename):
                return (self.lon == other.lon ) & (self.lat == other.lat)
            else:
                return False
        else:
            return False
    
    
    def read_raw(self, var):
        '''
        If profile is None the entire 2d array (depth,time) is returned
        If a profile object is provided, a 1d array (dimension depth) is returned
        corresponding to the pointprofile provided
        '''
        if not self.__file_already_read :
            ncIN = NC.netcdf_file(self.filename,'r')
            self._timesInFile = ncIN.variables['TIME'].data.copy()
            self._VAR         = ncIN.variables[var   ].data.copy()
            self._QC         = ncIN.variables[var + "_QC"  ].data.copy()
            if ncIN.variables.has_key('PRES'):
                self._PRES        = ncIN.variables['PRES'].data.copy()
                pres_qc           = ncIN.variables['PRES_QC'].data.copy()
            else:
                self._PRES        = ncIN.variables['DEPH'].data.copy()
                pres_qc           = ncIN.variables['DEPH_QC'].data.copy()

            self._TEMP = ncIN.variables['TEMP'].data.copy()
            self._PSAL = ncIN.variables['PSAL'].data.copy()

            temp_qc = ncIN.variables['TEMP_QC'].data.copy()
            psal_qc = ncIN.variables['PSAL_QC'].data.copy()
            ncIN.close()
            self._good_data = (temp_qc == 1 ) & (psal_qc == 1 ) & (pres_qc == 1 ) & (self._QC == 1)
            self.__file_already_read = True

            tconv = T90conv(self._TEMP)
            self._RHO = seawater.dens(self._PSAL,tconv, self._PRES)
        return self._VAR, self._PRES, self._QC, self._timesInFile

    def conversion(self,var,Profile,Rho):
        if var == 'DOX1':
            #from ml/l to mmol/m3
            return Profile*22.4 #grezza
        if var == 'CPHL':
            return Profile

    def read(self, var, time=None ):

        VAR, PRES, QC, timesInFile = self.read_raw(var)

        if time is None:
            return PRES,VAR  #returns array(times,depths)

        else:
            for it, t in enumerate(timesInFile):
                dateprofile = basetime + datetime.timedelta(days = t)
                if dateprofile == time:
                    Pres   = PRES[it,:]
                    Profile = VAR[it,:]
                    Qc      =  QC[it,:]
                    good    =  self._good_data[it,:]
                    rho     =  self._RHO[it,:]
                    break
            
            Pres    = Pres[good]
            Profile = Profile[good]
            rho     = rho[good]
            Qc      = Qc[good]
            Profile_inmodel_units = self.conversion(var, Profile, rho)
            return Pres, Profile_inmodel_units, Qc

    @staticmethod
    def fromfile(filename):

        nFiles=INDEX_FILE.size
        for iFile in range(nFiles):

            lon              = INDEX_FILE['geospatial_lon_min'][iFile]
            lat              = INDEX_FILE['geospatial_lat_min'][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile]
            available_params = INDEX_FILE['parameters'][iFile]
            if os.path.basename(filename) == os.path.basename(thefilename):
                return Mooring(lon,lat,filename,available_params)
        return None

def MooringSelector(var, T, region):

    LOC      = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/mooring/"

    INTEREST_VARLIST=['PHOS','SLCA','CPHL','NTRA','PHPH','DOX1']
    ALLVARLIST=['PHOS','SLCA','AMON','DOX1','DOX2','CPHL','NTRZ','NTRA','NTRI','PHPH']
    A = INDEX_FILE

    nFiles=A.size
    selected = []
    for i in range(nFiles):
        filename=os.path.basename(A['file_name'][i])
        parameters = A['parameters'][i]
        ti = TimeInterval(A[i]['time_coverage_start'],A[i]['time_coverage_end'],'%Y-%m-%dT%H:%M:%S' )
        lon = A[i]['geospatial_lon_min']
        lat = A[i]['geospatial_lat_min']

        if var is None:
            VarCondition = False
            for thevar in ALLVARLIST:
                if thevar in parameters: VarCondition = True
        else:
            VarCondition = var in parameters

        Filecondition = VarCondition & (filename[:2]=='MO') & ("mooring" in A['file_name'][i])
        
        if Filecondition & (ti.isInWindow(T)) & region.is_inside(lon, lat):
            filepath =  LOC + os.path.basename(A[i]['file_name'])
            ncIN = NC.netcdf_file(filepath,'r')
            timesInFile = ncIN.variables['TIME'].data.copy() # days since 1950-01-01T00:00:00Z
            ncIN.close()

            available_params = A['parameters'][i]
            M = Mooring(lon,lat,filepath,available_params)
            for t in timesInFile:
                dateprofile = basetime + datetime.timedelta(days = t)
                if T.contains(dateprofile):
                    mp = MooringProfile(dateprofile,lon,lat,M,available_params)
                    selected.append(mp)
    return selected
