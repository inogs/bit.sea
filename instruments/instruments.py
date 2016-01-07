
import bio_float
import mooring

from var_conversions import FLOATVARS, MOORINGVARS

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
    else:
        floatvar = FLOATVARS[var]
        mooringvar = MOORINGVARS[var]
    
    LIST=[]
    LIST.extend(bio_float.FloatSelector(floatvar  , T, region))
    LIST.extend(mooring.MooringSelector(mooringvar, T, region))
    return LIST

