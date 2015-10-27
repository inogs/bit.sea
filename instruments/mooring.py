import numpy as np
import os,datetime
from instruments.instrument import Instrument, Profile
import scipy.io.netcdf as NC

from commons.time_interval import TimeInterval
basetime = datetime.datetime(1950,1,1,0,0,0)

class MooringProfile(Profile):
    def __init__(self,time,lat,lon, my_mooring):
        self.time = time
        self.lat  = lat
        self.lon  = lon
        self.my_mooring = my_mooring

    def read(self,var):
        return self.my_mooring.read(var,self.time)

class Mooring(Instrument):
    def __init__(self, lon, lat, filename, params):
        self.lon = lat
        self.lat = lat
        self.filename = filename
        self.available_params = params
    
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
    
    def read(self, var, time=None ):

        ncIN = NC.netcdf_file(self.filename,'r')
        timesInFile = ncIN.variables['TIME'].data.copy()
        VAR         = ncIN.variables[var   ].data.copy()
        PRES        = ncIN.variables['PRES'].data.copy()
        ncIN.close()
        if time is None:
            return PRES,VAR  #returns array(times,depths)

        else:
            for it, t in enumerate(timesInFile):
                dateprofile = basetime + datetime.timedelta(days = t)
                if dateprofile == time:
                    Pres   = PRES[it,:]
                    Profile = VAR[it,:]
            return Pres,Profile

def MooringSelector(var, T, region):
    filename = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/index_monthly.txt"
    LOC      = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/COPERNICUS/mooring/"
    mydtype= np.dtype([('catalog_id','S20'),
              ('file_name','S200'),
              ('geospatial_lat_min',np.float32),
              ('geospatial_lat_max',np.float32),
              ('geospatial_lon_min',np.float32),
              ('geospatial_lon_max',np.float32),
              ('time_coverage_start','S19'),
              ('time_coverage_end','S19'),
              ('provider','S30'),
              ('date_update','S30'),
              ('data_mode','S1'),
              ('parameters','S200')] )

    ALLVARLIST=['PHOS','SLCA','AMON','DOX1','DOX2','CPHL','NTRZ','NTRA','NTRI','PHPH']
    INTEREST_VARLIST=['PHOS','SLCA','CPHL','NTRA','PHPH']
    #NTRZ:long_name = "nitrate + nitrite" ;

    A=np.loadtxt(filename,dtype=mydtype, delimiter=",", skiprows=6 )
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
            for thevar in INTEREST_VARLIST:
                if thevar in parameters: VarCondition = True
        else:
            VarCondition = var in parameters

        condition = VarCondition & (filename[:2]=='MO') & (ti.isInWindow(T)) & region.is_inside(lon, lat)
        if condition :
            filepath =  LOC + os.path.basename(A[i]['file_name'])
            ncIN = NC.netcdf_file(filepath,'r')
            timesInFile = ncIN.variables['TIME'].data.copy() # days since 1950-01-01T00:00:00Z
            ncIN.close()

            available_params = A['parameters'][i]

            for t in timesInFile:
                dateprofile = basetime + datetime.timedelta(days = t)
                if T.contains(dateprofile):
                    M = Mooring(lon,lat,filepath,available_params)
                    mp = MooringProfile(dateprofile,lon,lat,M)
                    selected.append(mp)
    return selected
