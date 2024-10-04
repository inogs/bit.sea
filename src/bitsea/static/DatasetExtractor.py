import numpy as np
import netCDF4 as NC4
import datetime
from bitsea.instruments.instrument import ContainerProfile
import seawater
from seawater.library import T90conv
from bitsea.commons.utils import find_index
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons import timerequestors

class DatasetExtractor():
    
    def __init__(self,filename, datasetname):
        with NC4.Dataset(filename,'r') as ncIN:
            self.DATA     = ncIN.variables['DATA'][:]
            self.UNITS    = [unit for unit in ncIN.variables['UNITS'][:]]
            self.VARIABLES= [variable for variable in ncIN.variables['VARIABLES'][:]]
            self.CRUISES  = [cruise for cruise in ncIN.variables['Cruises'][:]]
        self.datasetname = datasetname

    def unique_rows(self,data, prec=5):
        
        d_r = np.fix(data * 10 ** prec) / 10 ** prec + 0.0
        b = np.ascontiguousarray(d_r).view(np.dtype((np.void, d_r.dtype.itemsize * d_r.shape[1])))
        _, ia,ic = np.unique(b, return_index=True, return_inverse=True)
        return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ia, ic

    def profileGenerator(self,TIME,lon,lat,values,depth,dataset):
        '''
        Returns a profile list by reading an expanded table
        of TIME,lon,lat,values,depth,dataset
        
        '''
        nValues=len(TIME)
        M=np.zeros((nValues,3),dtype=np.float32)
        M[:,0] = TIME
        M[:,1] = lon
        M[:,2] = lat
        uniques,_,ib = self.unique_rows(M)
        
        nProfiles = len(uniques)
    
        Profilelist=[]
        for i in range(nProfiles):
            time = datetime.datetime.fromordinal(int(uniques[i,0]))
            lon  = uniques[i,1].astype(np.float32)
            lat  = uniques[i,2].astype(np.float32)
            inthisprofile = ib==i
            Values=values[inthisprofile]
            Depth = depth[inthisprofile]
            
            iddataset = int(dataset[inthisprofile][0])
            Cruise    = self.CRUISES[iddataset-1]
            LP = ContainerProfile(lon,lat,time,Depth,Values,Cruise, self.datasetname)
            Profilelist.append(LP)
        return Profilelist

    def selector(self,var,T_int, region):
        '''
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval or timerequestors.Clim_season object)
          region  (region object)

         For programmers : calculation of density should be moved in class constructor.
         '''
        assert isinstance(T_int, (TimeInterval, timerequestors.Clim_season, timerequestors.Clim_month))
        ivar  = find_index(var, self.VARIABLES)
        values= self.DATA[ivar,:].copy()
        units = self.UNITS[ivar]

        if False: #units =="\\mumol/kg":
            itemp  = find_index('temp'    , self.VARIABLES)
            ipsal  = find_index('salinity', self.VARIABLES)
            temp   = self.DATA[itemp,:]
            sali   = self.DATA[ipsal,:]
            pres   = self.DATA[5,:]
            good_rho  = (sali < 1.e+19 ) & (sali>0) & (temp < 1.e+19 ) & (temp>0) & (pres < 1.e+19 ) & (pres>0)
            t = T90conv(temp)
            n = len(values)
            calculated_rho  = np.ones((n),np.float32)*np.nan
            assumed_density = np.ones((n),np.float32)*np.nan
            calculated_rho[good_rho]  = seawater.dens(sali[good_rho],t[good_rho],pres[good_rho])

            for i in range(n):
                if good_rho[i]:
                    assumed_density[i] = calculated_rho[i]

            good = ~np.isnan(assumed_density)
            values[good] = values[good] * assumed_density[good] /1000.
            values[~good] = 1.e+20

        good  = (values < 1e+19) & (values > 0)
        
        values = values[good]
        
        year    = self.DATA[ 0,good].astype(int)
        month   = self.DATA[ 1,good].astype(int)
        day     = self.DATA[ 2,good].astype(int)
        lat     = self.DATA[ 3,good]
        lon     = self.DATA[ 4,good]
        depth   = self.DATA[ 5,good]
        dataset = self.DATA[-1,good]
        
        nValues=values.size
        Selected = np.zeros((nValues,),dtype=bool)
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
        
        if var == "pCO2_rec":
            Ptot= 1 + depth/10  # approximation for total pressure in atmosphere:
                                #! press atm + press water column (in atmosphere)
            values= values /  ( np.exp( ( 1-Ptot) *0.001366 ) )

        return self.profileGenerator(Time, Lon, Lat, values, depth, dataset)
    
    def cruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

         Returns a profile list
         
         '''
        iCruise = find_index(Cruisename, self.CRUISES)
        ivar    = find_index(var, self.VARIABLES)
        values= self.DATA[ivar,:]
        good = (values < 1e+19) & (values > 0) 
        ii = self.DATA[-1,:] == (iCruise+1)
        Selected = good & ii
        year    = self.DATA[ 0,Selected].astype(int)
        month   = self.DATA[ 1,Selected].astype(int)
        day     = self.DATA[ 2,Selected].astype(int)
        lat     = self.DATA[ 3,Selected]
        lon     = self.DATA[ 4,Selected]
        depth   = self.DATA[ 5,Selected]
        dataset = self.DATA[-1,Selected]
        
        values  = values[Selected]
        nValues=values.size
        TIME     = np.zeros((nValues), dtype=np.int32)
        for i in range(nValues):
            time = datetime.datetime(year[i],month[i],day[i])
            TIME[i] = time.toordinal()
        #print(iCruise,ivar,dataset,year)
        return self.profileGenerator(TIME, lon, lat, values, depth, dataset)
