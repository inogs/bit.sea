
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from DatasetExtractor import DatasetExtractor
from commons.utils import find_index

class CarbonReader(DatasetExtractor):
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Carbon/Dataset_Med_CarbSys_RIC.nc"
        self.DataExtractor = DatasetExtractor(self.filename,'Carbon')

        # DATA ELIMINATION in order to not duplicate values with nutrients dataset
        for cruisename in ['METEOR','METEOR51', 'METEOR95','PROSOPE']:
            iCruise = find_index(cruisename, self.DataExtractor.CRUISES)
            for var in ['nitrate','phosphate','silicate','oxygen']:
                ivar    = find_index(var, self.DataExtractor.VARIABLES)
                ii = self.DataExtractor.DATA[-1,:] == (iCruise+1)
                self.DataExtractor.DATA[ivar,ii] = -999.0

 
    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        var can be one of these:
         - DIC
         - ALK
         - temp
         - theta
         - salinity
         - silicate
         - nitrate
         - phosphate
         - oxygen
         - sigma_theta
         - density
         - sigma_t
         - pH25_sws
         - pH25_T
         - pCO2
         - PHt_{T-Press-ins}
         - PHsws_{Tins-0dB}
         - DICric
         - xCO2

         Cruisename can be one of these
            METEOR   METEOR51   METEOR95
            BOUM 2008
            PROSOPE
            EGEO APRIL          EGEO SEPT
            REGINA MARIS
            Garcia del Cid
            SESAME - ADRIATICO 2008
            CARBOGIB 1     CARBOGIB 2    CARBOGIB 3
            CARBOGIB 4     CARBOGIB 5    CARBOGIB 6
            GIFT 1      GIFT 2
            DYFAMED


         Returns a profile list
         
         '''
        return self.DataExtractor.cruiseSelector(var, Cruisename)
    
    def Selector(self,var,T_int, region):
        '''
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        var can be one of these:
         - DIC
         - ALK
         - temp
         - theta
         - salinity
         - silicate
         - nitrate
         - phosphate
         - oxygen
         - sigma_theta
         - density
         - sigma_t
         - pH25_sws
         - pH25_T
         - pCO2
         - PHt_{T-Press-ins}
         - PHsws_{Tins-0dB}
         - DICric
         - xCO2
         '''

        if var is None:
            Profilelist = list()
            for myvar in ['nitrate','phosphate','oxygen','silicate','DIC','ALK','temp','salinity','pH25_T','xCO2','pCO2']:
                Profilelist.extend(self.DataExtractor.selector(myvar, T_int, region))
            return Profilelist
        return self.DataExtractor.selector(var, T_int, region)
        



if __name__ == '__main__':
    
    var= 'nitrate';
    TI = TimeInterval('19900101','2005101','%Y%m%d')
    Reg= Rectangle(-6,36,30,46)
    C = CarbonReader()
    ProfileLIST = C.Selector('xCO2', TI, Reg)
    p = ProfileLIST[0]
    a = p.read('nitrate')
    print a
    
    ProfileLIST2 = C.CruiseSelector('nitrate', 'PROSOPE')
    LIST3 = C.Selector(None, TI, Reg)


        

