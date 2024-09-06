
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from DatasetExtractor import DatasetExtractor


class MassimiliReader():
    
    
    
    def __init__(self):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/MASSIMILI/Dataset_Massimili.nc"
        self.DataExtractor = DatasetExtractor(self.filename,'Massimili')


    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        var can be one of these:
         - CHL
         - nitrate
         - phosphate
         - silicate
         - oxygen

         Cruisename can be one of these
            'BOUSSOLE 2015'    'MOOSE'    'DYFAMED'    'BAM'

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
         - CHL
         - oxygen
         if var is None, no selection is done about variable
         '''
        if var is None:
            Profilelist=list()
            for myvar in ['nitrate','phosphate','CHL','oxygen']:
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
    TI = TimeInterval('2015','2016','%Y')
    Reg= Rectangle(-6,36,30,46)
    N = MassimiliReader()
    ProfileLIST = N.Selector('nitrate', TI, Reg)
    s=0
    for p in ProfileLIST:
        Pres,Profile,Qc = p.read('nitrate')
        s = s + len(Pres)
    print s
    
    #Cruisename='MOOSE'
    #ProfileLIST2 = N.CruiseSelector(var, Cruisename)


        

