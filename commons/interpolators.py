from scipy.interpolate import griddata
import numpy as np


def surf_interp_2d(Mask1, Mask2, Map2d):
    '''
    Interpolates a 2d numpy array over a new regular mesh
    Arguments:
    * Mask1, Mask2 * mask objects, defined
    * Map2d        * 2d numpy array defined on Mask1 
    
    Returns:
    * NearestMap * 2d numpy array defined on Mask2 
    Works correctly only in surface because it uses griddata
    
    '''
    tmask1 = Mask1.mask_at_level(0)
    tmask2 = Mask2.mask_at_level(0)
    nP = tmask1.sum()
    points = np.zeros((nP,2), np.float32)
     
    points[:,0] = Mask16.xlevels[tmask1]
    points[:,1] = Mask16.ylevels[tmask1]
    values = Map2d[tmask1]
    MAP2d_nearest =griddata( points, values, (Mask2.xlevels, Mask2.ylevels), 'nearest', fill_value=np.nan)
    MAP2d_nearest[~tmask2]=np.nan
    return MAP2d_nearest


if __name__ == "__main__":
    from commons import netcdf3    
    from commons.mask import Mask
    from commons.dataextractor import DataExtractor

    Mask24=Mask("/gpfs/scratch/userexternal/plazzari/eas_v12/eas_v12_6/wrkdir/MODEL/meshmask.nc",dzvarname="e3t_0")    
    Mask16=Mask('/gpfs/scratch/userexternal/gbolzon0/RA_COAST_03/wrkdir/MODEL/meshmask.nc')
    
    filename="/gpfs/scratch/userexternal/gbolzon0/RA_COAST_03/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/ave.20120816-12:00:00.P_l.nc"


    VAR = DataExtractor(Mask16,filename,"P_l").values
    VAR24 = surf_interp_2d(Mask16, Mask24, VAR[0,:,:])
    
    netcdf3.write_2d_file(VAR24, 'P_l', 'P_l_24.nc', Mask24, fillValue=1e+20)



