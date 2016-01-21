
import bio_float
import mooring
from static.Carbon_reader import CarbonReader
from static.Nutrients_reader import NutrientsReader

from var_conversions import FLOATVARS, MOORINGVARS, CARBONVARS, NUTRVARS

def Selector(var,T,region):
    '''
    Arguments:
       var is a string indicating variable,
          if var is None, no selection is done about variable
       T is as TimeInterval istance
       region is a region istance

    Returns
       a list of profile objects:
       they can be either BioFloatProfile or MooringProfile instances
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
    LIST.extend(mooring.MooringSelector(mooringvar, T, region))

    N = NutrientsReader()
    C = CarbonReader()
    LIST.extend(N.Selector(nutr_var, T, region))
    LIST.extend(C.Selector(carb_var, T, region))
    return LIST

