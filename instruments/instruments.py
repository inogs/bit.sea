
import bio_float
import mooring

FLOATVARS={'O2o':'DOXY', \
           'N3n':'NITRATE',  \
           'P_i':'CHLA'}   

MOORINGVARS={'O2o':'DOX1', \
             'N3n':'NOTFOUND',  \
             'P_i':'CPHL'}  

def Selector(var,T,region):
    
 
    
    LIST=[]
    LIST.append(bio_float.FloatSelector(  FLOATVARS[var], T, region))
    LIST.append(mooring.MooringSelector(MOORINGVARS[var], T, region))
    return LIST

