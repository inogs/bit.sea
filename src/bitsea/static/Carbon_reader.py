from os import PathLike
from pathlib import Path
from typing import Union

from bitsea.commons.time_interval import TimeInterval
from bitsea.basins.region import Rectangle
from bitsea.static.DatasetExtractor import DatasetExtractor
from bitsea.commons.utils import find_index, find_index_s


DEFAULT_FILENAME = Path(
    "/g100_scratch/userexternal/vdibiagi/EMODnet_2022/NEW_int/fromSC/Dataset_Med_Carbon_2023_NetCDF4.nc"
)


class CarbonReader(DatasetExtractor):
    
    DATA_VARS = (
        'nitrate', 'phosphate', 'oxygen', 'silicate', 'DIC', 'ALK', 'temp',
        'salinity', 'pH25_T', 'xCO2', 'pCO2', 'PHt_{T-Press-ins}_obs',
        'PHt_{T-Press-ins}_ric'
    )

    def __init__(self, filename: Union[str, PathLike] = DEFAULT_FILENAME):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename = filename
        self.DataExtractor = DatasetExtractor(self.filename,'Carbon')

        # DATA ELIMINATION in order to not duplicate values with nutrients dataset
        dataset = self.DataExtractor.DATA[-1,:]
        for cruisename in ['METEOR','METEOR51', 'METEOR95','PROSOPE','TALPRO','MEDWAVES','MSM72','GIANI']: 
            iCruise,namesC= find_index_s(cruisename,self.DataExtractor.CRUISES)
            for i in iCruise:
                ii = dataset==(i+1)
                for var in ['nitrate','phosphate','silicate','oxygen']:
                   ivar    = find_index(var, self.DataExtractor.VARIABLES)
                   self.DataExtractor.DATA[ivar,ii] = -999.0

        # deleting EMODNET values, identified by the variable flag
        iflag = find_index('flag', self.DataExtractor.VARIABLES)
        for var in ['phosphate','silicate','nitrate','oxygen']:
                ivar    = find_index(var, self.DataExtractor.VARIABLES)
                jj = self.DataExtractor.DATA[iflag,:] == 1
                self.DataExtractor.DATA[ivar,jj] = -999.0

    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list  by selecting for
        variable (string) and
        Cruisename (string)

        var can be one of these:
         - DIC
         - ALK
         - temp
         - salinity
         - silicate
         - nitrate
         - phosphate
         - oxygen
         - pressure
         - density
         - pH25_T
         - pCO2
         - DICric
         - xCO2
         - PHt_{T-Press-ins}_obs
         - PHt_{T-Press-ins}_ric

         Returns a profile list
         
         '''
        return self.DataExtractor.cruiseSelector(var, Cruisename)
    
    def Selector(self, var, T_int, region):
        """
        Returns a profile list by selecting for
          variable (string),
          T_int   (TimeInterval object)
          region  (region object)

        var can be one of these:
         - DIC
         - ALK
         - temp
         - salinity
         - silicate
         - nitrate
         - phosphate
         - oxygen
         - pressure
         - density
         - pH25_T
         - pCO2
         - DICric
         - xCO2
         - PHt_{T-Press-ins}_obs
         - PHt_{T-Press-ins}_ric
         """

        if var is None:
            profile_list = list()
            for myvar in self.DATA_VARS:
                profile_list.extend(
                    self.DataExtractor.selector(myvar, T_int, region)
                )
            return profile_list
        return self.DataExtractor.selector(var, T_int, region)
        



if __name__ == '__main__':
    
    var= 'nitrate'
    TI = TimeInterval('19990101','2005101','%Y%m%d')
    Reg= Rectangle(-6,36,30,46)
    C = CarbonReader()
    ProfileLIST = C.Selector('xCO2', TI, Reg)
    p = ProfileLIST[0]
    a = p.read('nitrate')
    print(a)
    
    ProfileLIST2 = C.CruiseSelector('nitrate', 'PROSOPE')
    LIST3 = C.Selector(None, TI, Reg)


        

