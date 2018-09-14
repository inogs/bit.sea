
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from DatasetExtractor import DatasetExtractor


class Maredat_reader():
    
    
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Maredat/Maredat.nc"
        self.DataExtractor = DatasetExtractor(self.filename)


    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        var can be one of these:
         - Chl_a
         - tot_Chl_a


         Cruisename can be one of these
         '001', '002', '003', ...'409'


         Returns a profile list
         
         '''
        return self.DataExtractor.cruiseSelector(var, Cruisename)
    
    def Selector(self,var,T_int, region):
        '''
        Data period: 1991 - 2005
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        can be one of these:
         - Chl_a
         - tot_Chl_a


         if var is None, no selection is done about variable
         '''
        if var is None:
            Profilelist=list()
            for myvar in ['Chl_a','tot_Chl_a']:
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
    
    var= 'Chl_a';
    TI = TimeInterval('2000','2001','%Y')
    Reg= Rectangle(0,20,30,46)
    N = Maredat_reader()
    ProfileLIST = N.Selector('Chl_a', TI, Reg)
    
    
    Cruisename='003'
    ProfileLIST2 = N.CruiseSelector(var, Cruisename)


        

