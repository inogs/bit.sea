from commons.time_interval import TimeInterval
from Nutrients_reader import NutrientsReader
from Carbon_reader import CarbonReader
from instruments.var_conversions import NUTRVARS
import numpy as np

N=NutrientsReader()
C=CarbonReader()



def DatasetInfo(modelvarname):
    '''
    Argument:
    * modelvarname * string, like 'N1p'
    Returns:
    * var     * variable name to access dataset
    * dataset * NutrientsReader or Carbonreader object
    '''
    if modelvarname in ['N1p','N3n','O2o','N5s']:
        var =  NUTRVARS[modelvarname]
        dataset     = N
    if modelvarname in ['O3h', 'Ac'] :
        var = 'ALK'
        dataset = C
    if modelvarname in ['O3c', 'DIC'] :
        var ='DIC'
        dataset = C
    if modelvarname not in ['N1p','N3n','O2o','N5s','O3h', 'Ac','O3c', 'DIC' ]:
        raise ValueError("variable not in static dataset ")
    return var, dataset



TI = TimeInterval("1950","2050","%Y")

def get_climatology(modelvarname, subbasinlist, LayerList ):
    '''
    Returns 
    * CLIM * [nsub, nLayers] numpy array, a basic annual climatology
    * STD  * [nsub, nLayers] numpy array, relative std values
    '''
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    CLIM    = np.zeros((nSub, nLayers), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers), np.float32)*np.nan
    var, Dataset = DatasetInfo(modelvarname)
    for isub, sub in enumerate(subbasinlist):
        Profilelist =Dataset.Selector(var, TI, sub)
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        for p in Profilelist: 
            pres, profile, _ = p.read(var)
            Pres   = np.concatenate((Pres,pres))
            Values = np.concatenate((Values,profile))
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<=layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean()
                STD[isub, ilayer] = Values[ii].std()
    return CLIM, STD