import numpy as np
import scipy.io.netcdf as NC
import datetime
from commons.time_interval import TimeInterval
from basins.region import Rectangle


class LeslieProfile():
    def __init__(self,lon,lat,time, depth, values,name):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.pres = depth
        self.profile= values
        self.name  = name

    def read(self,var):
        return self.pres, self.profile



def unique_rows(data, prec=5):
    
    d_r = np.fix(data * 10 ** prec) / 10 ** prec + 0.0
    b = np.ascontiguousarray(d_r).view(np.dtype((np.void, d_r.dtype.itemsize * d_r.shape[1])))
    _, ia,ic = np.unique(b, return_index=True, return_inverse=True)
    return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ia, ic


class NutrientsReader():
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Nutrients/Dataset_Med_Nutrients.nc"
        ncIN=NC.netcdf_file(filename,'r')
        self.DATA     = ncIN.variables['DATA'].data.copy()
        self.UNITS    = ncIN.variables['UNITS'].data.copy()
        self.VARIABLES= ncIN.variables['VARIABLES'].data.copy()
        self.CRUISES  = ncIN.variables['Cruises'].data.copy()
        ncIN.close()

    def _find_index(self, thestring, STRINGLIST):
        nStrings = STRINGLIST.shape[0]
        for istring in range(nStrings):
            strippedstring=STRINGLIST[istring,:].tostring().strip()
            if strippedstring == thestring: break
        else:
            print thestring + " Not Found"
            raise NameError
        return istring

    def _profileGenerator(self,TIME,lon,lat,values,depth,dataset):
        '''
        Returns a profile list by reading an expanded table
        of TIME,lon,lat,values,depth,dataset
        
        '''
        nValues=len(TIME)
        M=np.zeros((nValues,3),dtype=np.float32)
        M[:,0] = TIME
        M[:,1] = lon
        M[:,2] = lat
        uniques,_,ib = unique_rows(M)
        
        nProfiles = len(uniques)
    
        Profilelist=[]
        for i in range(nProfiles):
            time = datetime.datetime.fromordinal(int(uniques[i,0]))
            lon  = uniques[i,1].astype(np.float32)
            lat  = uniques[i,2].astype(np.float32)
            inthisprofile = ib==i
            Values=values[inthisprofile]
            Depth = depth[inthisprofile]
            
            iddataset = dataset[inthisprofile][0]
            Cruise    = self.CRUISES[iddataset-1,:].tostring().strip()
            LP = LeslieProfile(lon,lat,time,Depth,Values,Cruise)
            Profilelist.append(LP)
        return Profilelist



    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        var can be one of these:
         - nitrate
         - phosphate
         - silicate
         - oxygen

         Cruisename can be one of these
            06MT51/2  BIOPT06
            CANARI
            DYFAMED
            DYFAMED/PAPADOC - 99
            MEDCIESM
            MEDGOOS2  MEDGOOS3  MEDGOOS4  MEDGOOS5
            MELISSA 2004
            MT84_3
            NORBAL  NORBAL2 NORBAL3  NORBAL4
            POSEIDONE1M3A
            PROSOPE
            RHOFI 1   RHOFI 2   RHOFI 3
            SINAPSI-3   SINAPSI-4

         Returns a profile list
         
         '''
        iCruise = self._find_index(Cruisename, self.CRUISES)
        ivar    = self._find_index(var, self.VARIABLES)
        values= self.DATA[ivar,:]
        good = (values < 1e+19) & (values > 0) 
        ii = self.DATA[10,:] == (iCruise+1)
        Selected = good & ii
        year    = self.DATA[ 0,Selected]
        month   = self.DATA[ 1,Selected]
        day     = self.DATA[ 2,Selected]
        lat     = self.DATA[ 3,Selected]
        lon     = self.DATA[ 4,Selected]
        depth   = self.DATA[ 5,Selected]
        dataset = self.DATA[10,Selected]
        
        values  = values[Selected]
        nValues=values.size
        TIME     = np.zeros((nValues), dtype=np.int32)
        for i in range(nValues):
            time = datetime.datetime(year[i],month[i],day[i])
            TIME[i] = time.toordinal()
        
        return self._profileGenerator(TIME, lon, lat, values, depth, dataset)
    
    def Selector(self,var,T_int, region):
        '''
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        var can be one of these:
         - nitrate
         - phosphate
         - silicate
         - oxygen
         '''
        
        ivar  = self._find_index(var, self.VARIABLES)
        values= self.DATA[ivar,:]
        good  = (values < 1e+19) & (values > 0)
        
        values = values[good]
        
        year    = self.DATA[ 0,good]
        month   = self.DATA[ 1,good]
        day     = self.DATA[ 2,good]
        lat     = self.DATA[ 3,good]
        lon     = self.DATA[ 4,good]
        depth   = self.DATA[ 5,good]
        dataset = self.DATA[10,good]
        
        nValues=values.size
        
        Selected = np.zeros((nValues,),dtype=np.bool)
        TIME     = np.zeros((nValues), dtype=np.int32)
        for i in range(nValues):
            time = datetime.datetime(year[i],month[i],day[i])
            TIME[i] = time.toordinal()
            if T_int.contains(time) & region.is_inside(lon[i], lat[i]):
                Selected[i] = True 
    
        Time    = TIME[Selected]
        Lon     =  lon[Selected]
        Lat     =  lat[Selected]
        values  =  values[Selected]
        depth   =   depth[Selected]
        dataset = dataset[Selected]
        
        return self._profileGenerator(Time, Lon, Lat, values, depth, dataset)



if __name__ == '__main__':
    
    var= 'nitrate';
    TI = TimeInterval('20020101','20030101','%Y%m%d')
    Reg= Rectangle(0,20,30,46)
    A = NutrientsReader()
    ProfileLIST = A.Selector('nitrate', TI, Reg)
    
    
    Cruisename='MELISSA 2004'
    ProfileLIST2 = A.CruiseSelector(var, Cruisename)


        

