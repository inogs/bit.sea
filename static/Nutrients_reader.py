
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from DatasetExtractor import DatasetExtractor


class NutrientsReader():
    
    
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Nutrients/Dataset_Med_Nutrients.nc"
        self.DataExtractor = DatasetExtractor(self.filename)


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
        return self.DataExtractor.cruiseSelector(var, Cruisename)
    
    def Selector(self,var,T_int, region):
        '''
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        can be one of these:
         - nitrate
         - phosphate
         - silicate
         - oxygen
         if var is None, no selection is done about variable
         '''
        if var is None:
            Profilelist=list()
            for myvar in ['nitrate','phosphate','silicate','oxygen']:
                sublist=self.DataExtractor.selector(myvar, T_int, region)
                for p in sublist: 
                    if not p in Profilelist: Profilelist.append(p)
            return Profilelist
        else:
            return self.DataExtractor.selector(var, T_int, region)
        
    def getAllprofiles(self, T_INT):
        Reg = Rectangle(-6,36,30,46)
        return self.Selector(None,T_INT, Reg)
    

if __name__ == '__main__':
    
    var= 'nitrate';
    TI = TimeInterval('20020101','20030101','%Y%m%d')
    Reg= Rectangle(0,20,30,46)
    N = NutrientsReader()
    ProfileLIST = N.Selector('nitrate', TI, Reg)
    
    
    Cruisename='MELISSA 2004'
    ProfileLIST2 = N.CruiseSelector(var, Cruisename)


        

