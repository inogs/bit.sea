import numpy as np 
import seawater as sw 
from commons.dataextractor import DataExtractor

def get_density(filename, Maskobj):
    '''
    Arguments : 
    * filename * must contain vosaline and votemper, like T forcings or ave phys
    * Maskobj  * a Mask object, consistent with filename
    
    Returns: 
    * rho * numpy 3d array, consistent with Maskobj
    '''
 
    jpk, jpj, jpi = Maskobj.shape
    tmask = Maskobj.mask    
    PRES = np.zeros((jpk,jpj,jpi), np.float32)
    for k in range(jpk):
        PRES[k,:,:] = Maskobj.zlevels[k]
    
    
    TEMP = DataExtractor(Maskobj,filename,'votemper').values
    SALI = DataExtractor(Maskobj,filename,'vosaline').values
    
    RHO         = np.zeros((jpk,jpj,jpi), np.float32)
    T           = np.zeros((jpk,jpj,jpi), np.float32)
    T[tmask]    = sw.temp(SALI[tmask],TEMP[tmask],PRES[tmask])
    RHO[tmask]  = sw.dens(SALI[tmask],T[tmask],PRES[tmask])
    return RHO

if __name__=="__main__":
    from commons.mask import Mask
    maskfile="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc"
    TheMask = Mask(maskfile,dzvarname="e3t_0")
    filename="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/AVE_PHYS/ave.20150116-12:00:00.phys.nc"
    rho =get_density(filename, TheMask)
