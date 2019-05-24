
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from DatasetExtractor import DatasetExtractor


class Float_opt_reader():
    
    
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT/First_approach/Float_opt_dataset.nc"
        self.DataExtractor = DatasetExtractor(self.filename)


    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        var can be one of these:
         - chl
         - SAL
         - TEMP
         - Ed_380
         - Ed_412
         - Ed_490

         Cruisename can be one of these
          lovbio001i', 'lovbio015c', 'lovbio016c', 'lovbio016d', 'lovbio017b', 'lovbio018c', 'lovbio035b', 'lovbio039b', 'lovbio042c', 
          'lovbio053b', 'lovbio058c', 'lovbio063c', 'lovbio064b', 'lovbio064c', 'lovbio066c', 'lovbio066d', 'lovbio067c', 'lovbio068d', 
          'lovbio072c', 'lovbio083d', 'lovbio085d', 'lovbio088d', 'lovbio089d', 'lovbio090d', 'lovbio091d', 'lovbio093d', 'ogsbio001b', 
          'ogsbio002b', 'ogsbio003b', 'ogsbio004b', 'ogsbio006b'


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
         - chl
         - SAL
         - TEMP
         - Ed_380
         - Ed_412
         - Ed_490

         if var is None, no selection is done about variable
         '''
        if var is None:
            Profilelist=list()
            for myvar in ['chl','SAL','TEMP','Ed_490','Ed_380','Ed_412']:
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
    
    var= 'chl';
    TI = TimeInterval('20140101','20150101','%Y%m%d')
    Reg= Rectangle(0,20,30,46)
    N = Float_opt_reader()
    ProfileLIST = N.Selector('chl', TI, Reg)
    
    
    Cruisename='lovbio063c'
    ProfileLIST2 = N.CruiseSelector(var, Cruisename)


        

