from os import PathLike
from pathlib import Path
from typing import Union

from bitsea.commons.time_interval import TimeInterval
from bitsea.basins.region import Rectangle
from bitsea.static.DatasetExtractor import DatasetExtractor
from bitsea.commons.utils import find_index


DEFAULT_FILENAME = Path(
    "/g100_scratch/userexternal/vdibiagi/EMODnet_2022/NEW_int/fromSC/Dataset_Med_Nutrients_2023_NetCDF4.nc"
)


class NutrientsReader:

    DATA_VARS = (
        'nitrate', 'phosphate', 'silicate', 'oxygen', 'nitrite', 'ammonium',
        'chlorophyll', 'total_chlorophyll'
    )

    def __init__(self, filename: Union[str, PathLike] = DEFAULT_FILENAME):
        '''
        Reads the NetCDF Dataset
        '''
        self.filename = filename
        self.DataExtractor = DatasetExtractor(self.filename, 'Nutrients')

        # QC  section ----------------
        M = self.DataExtractor
        nvars, nData=M.DATA.shape
        selected = np.ones((nData,),bool)

        dataset = self.DataExtractor.DATA[-1,:]
        LIST_CHECK=['Barney','BIOPT06','Boussole']
        for inameCruise in LIST_CHECK:
          id_dataset,namesC= find_index_s(inameCruise,self.DataExtractor.CRUISES)
          for i in id_dataset:
            bad = dataset==(i+1)
            selected[bad] = False

        M.DATA = M.DATA[:,selected]

        self.DataExtractor.DATA = M.DATA




    def CruiseSelector(self, var,Cruisename):
        '''
        Returns a profile list by selecting for
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
            for myvar in self.DATA_VARS:
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
    from bitsea.basins import V2 as OGS
    import numpy as np
    var= 'nitrate'
    TI = TimeInterval('1998','2018','%Y')
    Reg= Rectangle(0,20,30,46)
    N = NutrientsReader()


    ProfileLIST = N.Selector('phosphate', TI, OGS.adr2)
    import matplotlib.pyplot as pl
    fig, ax = pl.subplots()

    for ip, p in enumerate(ProfileLIST):
        Pres, Value, Qc=p.read('phosphate')
        ax.plot(Value,Pres,'b.')

    ax.set_ylim(-10,1000)
    ax.set_xlabel('phos')
    ax.set_ylabel('depth')
    if not ax.yaxis_inverted():ax.invert_yaxis()
    fig.show()

    print(len(ProfileLIST))
    

    from bitsea.layer_integral import coastline
    c_lon,c_lat=coastline.get()
    Cruisename='BIOPT06'


    ProfileLIST2 = N.CruiseSelector('nitrate', Cruisename)
    nP = len(ProfileLIST2)
    Lon = np.zeros((nP), np.float32) * np.nan
    Lat = np.zeros((nP), np.float32) * np.nan
    for ip, p in enumerate(ProfileLIST2):
        Lon[ip]= p.lon
        Lat[ip] =p.lat
    fig,ax=pl.subplots()
    ax.plot(c_lon,c_lat, 'k')
    ax.plot(Lon,Lat,'b.')
    fig.show()


    fig, ax = pl.subplots()
    for p in ProfileLIST2:
        Pres, Value, Qc=p.read('phosphate')
        ax.plot(Value,Pres,'b.')
    ax.set_ylim(-10,1000)
    ax.set_xlabel('phos')
    ax.set_ylabel('depth')
    if not ax.yaxis_inverted(): ax.invert_yaxis()
    fig.show()

