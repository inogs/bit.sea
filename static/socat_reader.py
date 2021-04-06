
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from DatasetExtractor import DatasetExtractor
from commons.utils import find_index

class CO2_socat_reader(DatasetExtractor):
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/CO2_socat/SOCAT_INFO_FCO2.nc"
        #self.filename="/Users/gbolzon/Documents/workspace/PY/GP/SOCAT_INFO_FCO2.nc"
        self.DataExtractor = DatasetExtractor(self.filename,'Socat')


 
    def CruiseSelector(self, var,Cruisename):
        '''

         Returns a profile list
         '''
        return self.DataExtractor.cruiseSelector(var, Cruisename)

    def clim_month_selector(self,var,month,region):
        '''
        Specific case of Selector, but faster.
        Works only for climatologic month, between 2010 and 2016
        Returns: a list of values, depth is not taken in account

        '''
        ivar  = find_index(var, self.DataExtractor.VARIABLES)
        month_bool=self.DataExtractor.DATA[1,:]==month
        year=self.DataExtractor.DATA[0,:]
        year_bool = (year>= 2010) &  (year<=2016)
        ii = month_bool & year_bool
        data=self.DataExtractor.DATA[:,ii]
        _,nP=data.shape
        values=[]
        for i in range(nP):
            lat     = data[ 3,i]
            lon     = data[ 4,i]
            if region.is_inside(lon, lat):
                values.append(data[ivar,i])
        return values


    def Selector(self,var,T_int, region):
        '''
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        var can be one of these:
         -temp
         -sal
         fCO2
         -gvCO2        
         '''

        if var is None:
            Profilelist = list()
            for myvar in ['temp','sal','fCO2']:
                Profilelist.extend(self.DataExtractor.selector(myvar, T_int, region))
            return Profilelist
        return self.DataExtractor.selector(var, T_int, region)
        



if __name__ == '__main__':
    

    TI = TimeInterval('20150101','2018101','%Y%m%d')
    Reg= Rectangle(-6,36,30,46)
    C = CO2_socat_reader()
    ProfileLIST = C.Selector('fCO2', TI, Reg)
    p = ProfileLIST[0]
    a = p.read('fCO2')
    print a
    
    LIST3 = C.Selector(None, TI, Reg)


        

