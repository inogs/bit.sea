from bitsea.commons.time_interval import TimeInterval
from bitsea.commons import season
from bitsea.commons import timerequestors
from bitsea.static.Nutrients_reader import NutrientsReader
from bitsea.static.Carbon_reader import CarbonReader
from bitsea.instruments.var_conversions import NUTRVARS
import numpy as np
from bitsea.basins import V2 as OGS
from bitsea.basins.basin import SimplePolygonalBasin, ComposedBasin

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
    if modelvarname in ['N1p','N3n','O2o','N4n','N5s','P_l','N2n']:
        var =  NUTRVARS[modelvarname]
        dataset     = N
    elif modelvarname in ['O3h', 'Ac', 'ALK'] :
        var = 'ALK'
        dataset = C
    elif modelvarname in ['O3c', 'DIC'] :
        var ='DIC_merged'
        dataset = C
    elif modelvarname in ['pH', 'PH'] :
        var='pH_ins_merged'
        dataset = C
    elif modelvarname == 'pCO2' :
        var='pCO2_rec'
        dataset = C
    else:
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
    #if var in ["pH", "PH", "pCO2"] : return sub
    #
    #if (sub.name == "tyr1"):
    #    if var in ["O3c","O3h"]: return OGS.tyr2
    #    return OGS.tyr1

    #if (sub.name == "adr1"):
    #    if var in ["N1p","N3n","N5s"]: return OGS.adr2
    #    return OGS.adr1

    #if (sub.name == "ion1"):
    #    if var in ["N3n", "O3c", "O3h"]: return ComposedBasin('ion4', [OGS.swm2, OGS.ion2, OGS.tyr2], 'Neighbors of ion1')
    #    return OGS.ion1

    return sub

def QualityCheck(var,sub):
    '''
    returns True if expert evaluation tells to keep the data, False when it tells rejection
    '''
    if sub.name == 'adr1':
        return False

    if var == 'N4n': 
        if sub.name in ["alb","swm1","tyr1","tyr2","adr1","ion1","lev1","lev2","lev3","lev4"]:
            return False

    if var == 'O2o':
        if sub.name== 'tyr1':
            return False

    if var in ["O3c", "O3h", "DIC", "ALK", "pH", "PH", "PCO2"]:
        if sub.name in ["tyr1","adr1","lev3"]:
            return False

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

def mylogistic4(xdata,*x4):
    a = x4[0]
    b = x4[1]
    c = x4[2]
    d = x4[3]
    F = d + (a - d) / (1+ pow((xdata / c),b) )
    return F

# FIXME
TI = TimeInterval('19950101','20240101',"%Y%m%d")
#TI = TimeInterval('20060101','20240101',"%Y%m%d")

def get_climatology(modelvarname, subbasinlist, LayerList, useLogistic=False,startpt=np.asarray([0.1, 0.1, 500, 4],dtype=np.float64),basin_expand=False, QC=False,climseason=-1,climatology_interval=TI):
    '''
    basin_expand=True: apply research of data in close sub-basins (see basin_expansion)
    QC=True: exclude non-good subbasins (see QualityCheck)
    climseason=-1 (deafult): any season choice
    climseason=0,1,2,3: winter,spring,summer,autumn mean (standard seasons definition)
    useLogistic=True: (valid only for some nutrient data) it computes the climatology from
    a fit by a logistic curve on profiles only using data down to 1000m. Not recommended for climseason!=-1
    startpt must be an array of 4 values of np.float64
    useLogistic=False (default): it computes the climatology as means over layers
    Returns 
    * CLIM * [nsub, nLayers] numpy array, a basic or seasonal climatology, computed as mean or logistic on vertical layers
    * STD  * [nsub, nLayers] numpy array, relative std values
    * NVALS   * [nsub] numpy array, number of values in each layer
    * NPROFS  * [nsub] numpy array, number of profiles in each subbasin
    '''
    import scipy.optimize as opt
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    CLIM    = np.zeros((nSub, nLayers), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers), np.float32)*np.nan
    NVALS   = np.zeros((nSub, nLayers), np.int32)
    NPROFS  = np.zeros((nSub), np.int32)
    var, Dataset = DatasetInfo(modelvarname)
    var_exp      = Internal_conversion(modelvarname)
    p0_user=startpt

    if useLogistic==True:
        if p0_user.shape!=(4,):
            raise ValueError("the first guess for parameters is not an array of 4 values, as requested ")
        if p0_user.dtype!=np.float64:
            raise ValueError("the first guess for parameters is not an array of float64, as requested ")
        if modelvarname not in ['N1p','N3n','N5s']:
            raise ValueError("variable not properly interpolable by a logistic curve (as instead N1p, N3n, N5s)")
        depthLayer = np.zeros(nLayers, np.float32)*np.nan
        for ilayer, layer in enumerate(LayerList):
            depthLayer[ilayer] = 0.5*(layer.top+layer.bottom)

    if climseason==-1:
        T_int = climatology_interval
    if climseason in [0,1,2,3]:
        S=season.season()
        T_int = timerequestors.Clim_season(climseason,S)
    for isub, sub in enumerate(subbasinlist):
        ProfilelistAll =Dataset.Selector(var, T_int, sub)
        Profilelist=[p for p in ProfilelistAll if TI.contains(p.time)]
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        nProfiles=0
        for p in Profilelist: 
            pres, profile, _ = p.read(var)
            Pres   = np.concatenate((Pres,pres))
            Values = np.concatenate((Values,profile))
            nProfiles=nProfiles+1
        NPROFS[isub] = nProfiles    
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean() 
                STD[isub, ilayer] = Values[ii].std()
                NVALS[isub,ilayer] = ii.sum()
                # CLIM from logistic curve down to 1000 m
        if useLogistic:
            Pres_fit = Pres[Pres<=1000]
            Values_fit = Values[Pres<=1000]
            if nProfiles>0:
                try:
                    popt, pcov = opt.curve_fit(mylogistic4,Pres_fit,Values_fit,p0=p0_user,method="lm")
                    #print("popt: ",popt)
                    clim = mylogistic4(depthLayer, *popt)
                except RuntimeError:
                    print("Error - curve_fit failed")
                    clim=np.nan
            else:
                clim=np.nan
            CLIM[isub,:] = clim

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
            print(INDEX_LIST)
            CLIM[isub,:] = CLIM[INDEX_LIST,:].mean(axis=0)
            STD [isub,:] =  STD[INDEX_LIST,:].mean(axis=0) # brutto ...
            NVALS[isub,:] = 0
            NPROFS[isub] = 0
    if QC:
        for isub, sub in enumerate(subbasinlist):
            is_good = QualityCheck(var_exp, sub)
            if not is_good:
                CLIM[isub,:] = np.nan
                STD[isub,:] = np.nan
                NVALS[isub,:] = 0
                NPROFS[isub] = 0         

    return CLIM, STD, NVALS, NPROFS


def get_climatology_open(modelvarname, subbasinlist, LayerList, TheMask, limdepth=200.0, obsdepthlim=50.0,useLogistic=False,startpt=np.asarray([0.1, 0.1, 500, 4],dtype=np.float64),basin_expand=False, QC=False,climseason=-1, climatology_interval=TI):
    '''
    get climatology excluding coastal areas, it needs TheMask (a mask object)
    limit depth for coastal areas (limdepth) is 200 by default
    It consider only profiles with data deeper than obsdepthlim (default 50 m) 
    basin_expand=True: apply research of data in close sub-basins (see basin_expansion)
    QC=True: exclude non-good subbasins (see QualityCheck)
    climseason=-1 (default): any season choice
    climseason=0,1,2,3: winter,spring,summer,autumn mean (standard seasons definition)
    useLogistic=True: (valid only for some nutrient data) it computes the climatology from 
    a fit by a logistic curve on profiles only using data down to 1000m.Not recommended for climseason!=-1
    startpt must be an array of 4 values of np.float64
    useLogistic=False (default): it computes the climatology as means over layers
    Returns 
    * CLIM * [nsub, nLayers] numpy array, a basic or seasonal climatology, computed as mean or logistic on vertical layers
    * STD  * [nsub, nLayers] numpy array, always computed as std on vertical layers
    * NVALS   * [nsub] numpy array, number of values in each layer
    * NPROFS  * [nsub] numpy array, number of profiles in each subbasin
    '''

    import scipy.optimize as opt
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    CLIM    = np.zeros((nSub, nLayers), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers), np.float32)*np.nan
    NVALS   = np.zeros((nSub, nLayers), np.int32)
    NPROFS  = np.zeros((nSub), np.int32)
    var, Dataset = DatasetInfo(modelvarname)
    var_exp      = Internal_conversion(modelvarname)
    mask200 = TheMask.mask_at_level(limdepth)
    p0_user=startpt
    #print(p0_user)

    if useLogistic==True:
        if p0_user.shape!=(4,):
            raise ValueError("the first guess for parameters is not an array of 4 values, as requested ")
        if p0_user.dtype!=np.float64:
            raise ValueError("the first guess for parameters is not an array of float64, as requested ")
        if modelvarname not in ['N1p','N3n','N5s']:
            raise ValueError("variable not properly interpolable by a logistic curve (as instead N1p, N3n, N5s)")
        depthLayer = np.zeros(nLayers, np.float32)*np.nan
        for ilayer, layer in enumerate(LayerList):
            depthLayer[ilayer] = 0.5*(layer.top+layer.bottom)

    if climseason==-1:
        T_int = climatology_interval
    if climseason in [0,1,2,3]:
        S=season.season()
        T_int = timerequestors.Clim_season(climseason,S)
    for isub, sub in enumerate(subbasinlist):
        print(sub)
        ProfilelistAll =Dataset.Selector(var, T_int, sub) 
        Profilelist=[p for p in ProfilelistAll if TI.contains(p.time)]
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        nProfiles=0 # in the subd
        for p in Profilelist:
            lonp = p.lon
            latp = p.lat
            timep = p.time
            ilon,ilat = TheMask.convert_lon_lat_to_indices(lonp,latp)
            if mask200[ilat,ilon]:
                pres, profile, _ = p.read(var)
                if max(pres)>=obsdepthlim: # to delete observations that are only in the surface layers (not real profiles)
                    Pres   = np.concatenate((Pres,pres))
                    Values = np.concatenate((Values,profile))
                    nProfiles=nProfiles+1
        NPROFS[isub] = nProfiles
        # CLIM and STD as statistics on vertical layers
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean()
                STD[isub, ilayer] = Values[ii].std()
                NVALS[isub,ilayer] = ii.sum()
        # CLIM from logistic curve down to 1000 m
        if useLogistic:
            Pres_fit = Pres[Pres<=1000]
            Values_fit = Values[Pres<=1000]
            if nProfiles>0:
                try:
                    popt, pcov = opt.curve_fit(mylogistic4,Pres_fit,Values_fit,p0=p0_user,method="lm")
                    #print("popt: ",popt)
                    clim = mylogistic4(depthLayer, *popt)
                except RuntimeError:
                    print("Error - curve_fit failed")
                    clim=np.nan
            else:
                clim=np.nan
            CLIM[isub,:] = clim

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
            print(INDEX_LIST)
            CLIM[isub,:] = CLIM[INDEX_LIST,:].mean(axis=0)
            STD [isub,:] =  STD[INDEX_LIST,:].mean(axis=0) # brutto ...
            NVALS[isub,:] = 0
            NPROFS[isub] = 0

    if QC:
        for isub, sub in enumerate(subbasinlist):
            is_good = QualityCheck(var_exp, sub)
            if not is_good:
                CLIM[isub,:] = np.nan
                STD[isub,:] = np.nan
                NVALS[isub,:] = 0 
                NPROFS[isub] = 0
    
    return CLIM, STD, NVALS, NPROFS
    

def get_histo_open(modelvarname, subbasinlist, LayerList, TheMask, limdepth = 200.0, nbins=50, climatology_interval=TI):
    '''
    get histogram of values in open sea areas (number of values VS values)
    returns: intervals of values of the variable, number of values
    * VALUES * [nsub, nLayers, nbin+1]
    * NUMBER * [nsub, nLayers, nbin]
    '''
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    VALUES     = np.zeros((nSub, nLayers, nbins+1), np.float32)*np.nan
    NUMBER     = np.zeros((nSub, nLayers, nbins), np.float32)*np.nan
    var, Dataset = DatasetInfo(modelvarname)
    var_exp      = Internal_conversion(modelvarname)
    mask200 = TheMask.mask_at_level(limdepth)

    T_int = climatology_interval
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
            if ii.sum()> 0:
                TMP=np.histogram(Values[ii], nbins)
                NUMBER[isub,ilayer,:]=TMP[0]
                VALUES[isub,ilayer,:]=TMP[1]

    return NUMBER, VALUES

def get_climatology_coast(modelvarname, subbasinlist, LayerList, TheMask, limdepth=200.0, basin_expand=False, QC=False, climseason=-1, climatology_interval=TI):
    '''
    get climatolgy in coastal areas, it needs TheMask (a mask object)
    limit depth for coastal areas (limdepth) is 200 by default
    basin_expand=True: apply research of data in close sub-basins (see basin_expansion)
    QC=True: exclude non-good subbasins (see QualityCheck)
    climseason=-1 (deafult): any season choice
    climseason=0,1,2,3: winter,spring,summer,autumn mean (standard seasons definition)
    Returns 
    * CLIM * [nsub, nLayers] numpy array, a basic or seasonal climatology
    * STD  * [nsub, nLayers] numpy array, relative std values
    * NVALS   * [nsub] numpy array, number of values in each layer
    * NPROFS  * [nsub] numpy array, number of profiles in each subbasin
    '''
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    CLIM    = np.zeros((nSub, nLayers), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers), np.float32)*np.nan
    NVALS   = np.zeros((nSub, nLayers), np.int32)
    NPROFS  = np.zeros((nSub), np.int32)
    var, Dataset = DatasetInfo(modelvarname)
    var_exp      = Internal_conversion(modelvarname)
    masklimdepth = TheMask.mask_at_level(limdepth)
    mask0 = TheMask.mask_at_level(0.0)
    coastmask=mask0 & (~masklimdepth)

    if climseason==-1:
        T_int = climatology_interval
    if climseason in [0,1,2,3]:
        S=season.season()
        T_int = timerequestors.Clim_season(climseason,S)
    for isub, sub in enumerate(subbasinlist):
        ProfilelistAll =Dataset.Selector(var, T_int, sub)
        Profilelist=[p for p in ProfilelistAll if TI.contains(p.time)]
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        nProfiles=0 
        for p in Profilelist:
            lonp = p.lon
            latp = p.lat
            ilon,ilat = TheMask.convert_lon_lat_to_indices(lonp,latp)
            if coastmask[ilat,ilon]:
                pres, profile, _ = p.read(var)
                Pres   = np.concatenate((Pres,pres))
                Values = np.concatenate((Values,profile))
                nProfiles=nProfiles+1
        NPROFS[isub] = nProfiles        
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean()
                STD[isub, ilayer] = Values[ii].std()
                NVALS[isub,ilayer] = ii.sum()

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
            print(INDEX_LIST)
            CLIM[isub,:] = CLIM[INDEX_LIST,:].mean(axis=0)
            STD [isub,:] =  STD[INDEX_LIST,:].mean(axis=0) # brutto ...
            NVALS[isub,:] = 0
            NPROFS[isub] = 0

    if QC:
        for isub, sub in enumerate(subbasinlist):
            is_good = QualityCheck(var_exp, sub)
            if not is_good:
                CLIM[isub,:] = np.nan
                STD [isub,:] = np.nan
                NVALS[isub,:] = 0
                NPROFS[isub] = 0

    return CLIM, STD, NVALS, NPROFS


if __name__ == "__main__":
    from bitsea.commons.layer import Layer
    PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500,3000,4000,5000])
    LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]
    SUBLIST = OGS.P.basin_list
    N1p_clim, N1p_std, _, _ = get_climatology('N1p', SUBLIST, LayerList, basin_expand=False, QC=True)
    N5s_clim, N5s_std, _, _ = get_climatology('N5s', SUBLIST, LayerList, basin_expand=False, QC=True)
    N4n_clim, N4n_std, _, _ = get_climatology('N4n', SUBLIST, LayerList, basin_expand=False, QC=True)
