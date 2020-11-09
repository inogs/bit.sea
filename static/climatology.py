from commons.time_interval import TimeInterval
from commons import season
from commons import timerequestors
from Nutrients_reader import NutrientsReader
from Carbon_reader import CarbonReader
from instruments.var_conversions import NUTRVARS
import numpy as np
from basins import V2 as OGS
from basins.basin import SimplePolygonalBasin, ComposedBasin

N=NutrientsReader()
C=CarbonReader()



def DatasetInfo(modelvarname):
    '''
    1) Helps the user to navigate inside static dataset,
    because some variables is in Nutrients one, some other
    are in the Carbon one.
    2) Defines the official name of variables to be used in deliverables
       statistics.

    Argument:
    * modelvarname * string, like 'N1p'
    Returns:
    * var     * variable name to access dataset
    * dataset * NutrientsReader or Carbonreader object
    '''
    if modelvarname in ['N1p','N3n','O2o','N4n','N5s']:
        var =  NUTRVARS[modelvarname]
        dataset     = N
    if modelvarname in ['O3h', 'Ac', 'ALK'] :
        var = 'ALK'
        dataset = C
    if modelvarname in ['O3c', 'DIC'] :
        var ='DICric'
        dataset = C
    if modelvarname in ['pH', 'PH'] :
        var='PHt_{T-Press-ins}'
        dataset = C
    if modelvarname == 'pCO2' :
        var='pCO2'
        dataset = C
    if modelvarname not in ['N1p','N3n','O2o','N4n','N5s','O3h', 'Ac','ALK','O3c', 'DIC', 'pH',',PH', 'pCO2' ]:
        raise ValueError("variable not in static dataset ")
    return var, dataset

def Internal_conversion(modelvarname):
    if modelvarname in ['O3c', 'DIC'] : return "O3c"
    if modelvarname in ['O3h', 'Ac', 'ALK'] : return "O3h"
    return modelvarname


def basin_expansion(sub, var):
    '''
    In order to have more data for each subbasin, sometimes we need to
    expand the research to nearest neighbours.
    The amount of available data depends on subbasin and var.
    This function follows the table provided by Valeria Di Biagio
    and returns a suitable subbasin (based on nearest neighbours) to guarantee a
    minimum amount of data in research.

    Arguments:
    * sub * a basin object
    * var * string, can be "N1p","N3n","N5s","O2o","O3c","O3h"
    Returns:
    * search_sub * a sub object
    '''
    #assert var in ["N1p","N3n","N5s","O2o","O3c","O3h"]
    if var in ["pH", "PH", "pCO2"] : return sub

    if (sub.name == "tyr1"):
        if var in ["O3c","O3h"]: return OGS.tyr2
        return OGS.tyr1

    if (sub.name == "adr1"):
        if var in ["N1p","N3n","N5s"]: return OGS.adr2
        return OGS.adr1

    if (sub.name == "ion1"):
        if var in ["N3n", "O3c", "O3h"]: return ComposedBasin('ion4', [OGS.swm2, OGS.ion2, OGS.tyr2], 'Neighbors of ion1')
        return OGS.ion1

    return sub

def QualityCheck(var,sub):
    '''
    returns True if expert evaluation tells to keep the data, False when it tells rejection
    '''
    if sub.name == 'adr1':
        if var in ["O2o", "N3n"] : return False
    if sub.name =="tyr1":
        if var in ["O2o", "O3c", "O3h"]: return False
    return True

def get_sub_indexes(sub_search):
    '''
    Returns a list of indexes
    '''
    if isinstance(sub_search, SimplePolygonalBasin):
        index = OGS.Pred.basin_list.index(sub_search)
        return [index]
    if isinstance(sub_search, ComposedBasin):
        L=[]
        for sub in sub_search:
            L.append(OGS.Pred.basin_list.index(sub))
        return L



TI = TimeInterval("1997","2016","%Y")

def get_climatology(modelvarname, subbasinlist, LayerList, basin_expand=False, QC=False,climseason=-1):
    '''
    basin_expand=True: apply research of data in close sub-basins (see basin_expansion)
    QC=True: exclude non-good subbasins (see QualityCheck)
    climseason=-1 (deafult): any season choice
    climseason=0,1,2,3: winter,spring,summer,autumn mean (standard seasons definition)
    Returns 
    * CLIM * [nsub, nLayers] numpy array, a basic or seasonal climatology
    * STD  * [nsub, nLayers] numpy array, relative std values
    '''
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    CLIM    = np.zeros((nSub, nLayers), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers), np.float32)*np.nan
    var, Dataset = DatasetInfo(modelvarname)
    var_exp      = Internal_conversion(modelvarname)
    if climseason==-1:
        T_int = TI
    if climseason in [0,1,2,3]:
        S=season.season()
        T_int = timerequestors.Clim_season(climseason,S)
    for isub, sub in enumerate(subbasinlist):
        Profilelist =Dataset.Selector(var, T_int, sub)
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        for p in Profilelist: 
            pres, profile, _ = p.read(var)
            Pres   = np.concatenate((Pres,pres))
            Values = np.concatenate((Values,profile))
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean()
                STD[isub, ilayer] = Values[ii].std()

    if basin_expand:
# 1. elimination nans by interpolation
        nLayers=len(LayerList)
        Layer_center=np.zeros((nLayers,),np.float32)
        for i,l in enumerate(LayerList): Layer_center[i] = (l.top + l.bottom)/2
        for isub, sub in enumerate(subbasinlist):
            y = CLIM[isub,:]
            nans= np.isnan(y)
            z_good = Layer_center[~nans]
            y_good = y[~nans]
            if len(y_good) > 0:
                CLIM[isub,:] = np.interp(Layer_center, z_good, y_good).astype(np.float32)
# 2 apply expansion following Valeria's table
        for isub, sub in enumerate(subbasinlist):
            sub_search = basin_expansion(sub, var_exp)
            INDEX_LIST=get_sub_indexes(sub_search)
            print INDEX_LIST
            CLIM[isub,:] = CLIM[INDEX_LIST,:].mean(axis=0)
            STD [isub,:] =  STD[INDEX_LIST,:].mean(axis=0) # brutto ...
    if QC:
        for isub, sub in enumerate(subbasinlist):
            is_good = QualityCheck(var_exp, sub)
            if not is_good:
                CLIM[isub,:] = np.nan
                STD [isub,:] = np.nan

    return CLIM, STD


def get_climatology_open(modelvarname, subbasinlist, LayerList, TheMask, limdepth=200, basin_expand=False, QC=False,climseason=-1):
    '''
    get climatolgy excluding coastal areas, it needs TheMask (a mask object)
    limit depth for coastal areas (limdepth) is 200 by default
    basin_expand=True: apply research of data in close sub-basins (see basin_expansion)
    QC=True: exclude non-good subbasins (see QualityCheck)
    climseason=-1 (deafult): any season choice
    climseason=0,1,2,3: winter,spring,summer,autumn mean (standard seasons definition)
    Returns 
    * CLIM * [nsub, nLayers] numpy array, a basic or seasonal climatology
    * STD  * [nsub, nLayers] numpy array, relative std values
    '''
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    CLIM    = np.zeros((nSub, nLayers), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers), np.float32)*np.nan
    var, Dataset = DatasetInfo(modelvarname)
    var_exp      = Internal_conversion(modelvarname)
    mask200 = TheMask.mask_at_level(limdepth)

    if climseason==-1:
        T_int = TI
    if climseason in [0,1,2,3]:
        S=season.season()
        T_int = timerequestors.Clim_season(climseason,S)
    for isub, sub in enumerate(subbasinlist):
        Profilelist =Dataset.Selector(var, T_int, sub)
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        for p in Profilelist: 
            lonp = p.lon
            latp = p.lat
            ilon,ilat = TheMask.convert_lon_lat_to_indices(lonp,latp)
            if mask200[ilat,ilon]:
                pres, profile, _ = p.read(var)
                Pres   = np.concatenate((Pres,pres))
                Values = np.concatenate((Values,profile))
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean()
                STD[isub, ilayer] = Values[ii].std()

    if basin_expand:
# 1. elimination nans by interpolation
        nLayers=len(LayerList)
        Layer_center=np.zeros((nLayers,),np.float32)
        for i,l in enumerate(LayerList): Layer_center[i] = (l.top + l.bottom)/2
        for isub, sub in enumerate(subbasinlist):
            y = CLIM[isub,:]
            nans= np.isnan(y)
            z_good = Layer_center[~nans]
            y_good = y[~nans]
            if len(y_good) > 0:
                CLIM[isub,:] = np.interp(Layer_center, z_good, y_good).astype(np.float32)
# 2 apply expansion following Valeria's table
        for isub, sub in enumerate(subbasinlist):
            sub_search = basin_expansion(sub, var_exp)
            INDEX_LIST=get_sub_indexes(sub_search)
            print INDEX_LIST
            CLIM[isub,:] = CLIM[INDEX_LIST,:].mean(axis=0)
            STD [isub,:] =  STD[INDEX_LIST,:].mean(axis=0) # brutto ...
    if QC:
        for isub, sub in enumerate(subbasinlist):
            is_good = QualityCheck(var_exp, sub)
            if not is_good:
                CLIM[isub,:] = np.nan
                STD [isub,:] = np.nan

    return CLIM, STD




if __name__ == "__main__":
    from commons.layer import Layer
    PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500,3000,4000,5000])
    LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]
    SUBLIST = OGS.P.basin_list
    N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList, basin_expand=True, QC=True)
    N5s_clim, N5s_std = get_climatology('N5s', SUBLIST, LayerList, basin_expand=True, QC=True)
    N4n_clim, N4n_std = get_climatology('N4n', SUBLIST, LayerList, basin_expand=True, QC=True)
