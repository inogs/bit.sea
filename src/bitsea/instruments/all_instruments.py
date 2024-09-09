import bio_float
import lovbio_float
import superfloat
import optbio_float
import optbio_float_2019
#import mooring
from static.Carbon_reader import CarbonReader
from static.Nutrients_reader import NutrientsReader

from var_conversions import LOVFLOATVARS, FLOATVARS, MOORINGVARS, CARBONVARS, NUTRVARS, FLOAT_OPT_VARS, FLOAT_OPT_VARS_2019, SUPERFLOAT_VARS

def Selector(var,T,region):
    '''
    Arguments:
       var is a string indicating variable,
          if var is None, no selection is done about variable
       T is as TimeInterval istance
       region is a region istance

    Returns
       a list of profile objects:
       they can be
        - BioFloatProfile or
        - MooringProfile or
        - ContainerProfile instances
    '''

    if var is None:
        floatvar = None
        mooringvar = None
        nutr_var   = None
        carb_var   = None
    else:
        floatvar = FLOATVARS[var]
        mooringvar = MOORINGVARS[var]
        nutr_var   = NUTRVARS[var]
        carb_var   = CARBONVARS[var]

    LIST=[]
    LIST.extend(bio_float.FloatSelector(floatvar  , T, region))
    #LIST.extend(mooring.MooringSelector(mooringvar, T, region))

    N = NutrientsReader()
    C = CarbonReader()
    LIST.extend(N.Selector(nutr_var, T, region))
    LIST.extend(C.Selector(carb_var, T, region))
    return LIST

def static_Selector(var,T,region):
    '''
    Arguments:
       var is a string indicating variable,
          if var is None, no selection is done about variable
       T is as TimeInterval istance
       region is a region istance

    Returns
       a list of CaontainerProfile instances :

    '''

    if var is None:
        nutr_var   = None
        carb_var   = None
    else:
        nutr_var   = NUTRVARS[var]
        carb_var   = CARBONVARS[var]

    LIST=[]

    N = NutrientsReader()
    C = CarbonReader()
    LIST.extend(N.Selector(nutr_var, T, region))
    LIST.extend(C.Selector(carb_var, T, region))
    return LIST

def getAllProfiles(TI):
    from basins.region import Rectangle
    All_Med = Rectangle(-6,36,30,46)
    return Selector(None, TI, All_Med)
