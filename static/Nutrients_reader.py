
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
         - nitrite
         - ammonium
         - phosphate
         - silicate
         - oxygen
         - chlorophyll
         - total_chlorophyll

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
            AIRWIN BARMED Barney BEHEMOTH
            BIOPRHOFI BOUSSOLE CASCADE CHACCRA
            COSIMO15 CYBO DEEP DYFAMED ECOLOPHY
            EMTEC ESTIME EUROSITES FLIPER GEOTETHYS
            GOLTS GYROSCOP HaiSec HERMES HIVERN
            HYGAM05 INTERREG ISOFLORE JELLYWATCH
            JUVALION LATEX MATER METROMED MINERCOT
            MODELFOS MOLA MOOGLI MOOSE MTPII-MATER
            NICOP OMER OPERA PRIMI PROPECHE SESAME
            SESIL SOFI STRATA.PRODELTA UNIMED
            WB13 WB14 Yoyo



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
         - nitrite
         - ammonium
         - phosphate
         - silicate
         - oxygen
         - chlorophyll
         - total_chlorophyll
         if var is None, no selection is done about variable
         '''
        if var is None:
            Profilelist=list()
            for myvar in ['nitrate','phosphate','silicate','oxygen','nitrite','ammonium','chorophyll','total_chlorophyll']:
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
    from basins import V2 as OGS
    var= 'nitrate';
    TI = TimeInterval('1998','2018','%Y')
    Reg= Rectangle(0,20,30,46)
    N = NutrientsReader()

    ProfileLIST = N.Selector('nitrate', TI, OGS.aeg)
    print len(ProfileLIST)
    
    
    Cruisename='MELISSA 2004'
    ProfileLIST2 = N.CruiseSelector(var, Cruisename)


        

