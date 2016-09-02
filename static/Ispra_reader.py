
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from IspraExtractor import IspraExtractor


class IspraReader():
    
    
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
#        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Nutrients/Dataset_Med_Nutrients.nc"
        self.filename="/pico/scratch/userexternal/ateruzzi/ispra2013.nc"
        self.DataExtractor = IspraExtractor(self.filename)


    def StationSelector(self, var, stationname):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        For variable names see the list in the original file

        For station names see the original file

        Returns a profile list
         
         '''
        return self.DataExtractor.stationSelector(var, stationname)
    
    def Selector(self,var,T_int, region):
        '''
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        For variable names see the original file
         if var is None only chlorophyll is extracted
         '''
        if var is None:
            Profilelist=list()
            myvar = 'chlrophyll a'
            Profilelist.extend(self.DataExtractor.selector(myvar, T_int, region))
            return Profilelist
        else:
            return self.DataExtractor.selector(var, T_int, region)

if __name__ == '__main__':
    
    var= 'nitrate';
    TI = TimeInterval('20130101','20130108','%Y%m%d')
    Reg= Rectangle(0,20,30,45)
    N = IspraReader()
    ProfileLIST = N.Selector(var, TI, Reg)
    
    
#    stationname='IT19CW0815301'
#    ProfileLIST2 = N.StationSelector(var, stationname)


        

