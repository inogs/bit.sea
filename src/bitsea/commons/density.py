import numpy as np
import gsw
from bitsea.commons.dataextractor import DataExtractor
from pathlib import Path
from bitsea.timeseries.plot import read_pickle_file
from bitsea.commons.mask import Mask

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
    LON         = np.broadcast_to(Maskobj.xlevels, (jpk,jpj,jpi))
    LAT         = np.broadcast_to(Maskobj.ylevels, (jpk,jpj,jpi))
    SA          = gsw.SA_from_SP(SALI[tmask], PRES[tmask], LON[tmask], LAT[tmask])
    CT          = gsw.CT_from_pt(SA, TEMP[tmask])
    RHO[tmask]  = gsw.rho(SA, CT, PRES[tmask])
    return RHO


def get_density_from_4D_profiles(T, S, Maskobj: Mask):
    '''
    Arguments :
    * T4D * is a 4D (nFrames,nSub,nCoast,jpk) array of temperature
    * S4D * is a 4D (nFrames,nSub,nCoast,jpk) array of salinity
    * Maskobj  * a Mask object, consistent with filename

    Returns:
    * rho * numpy 3d array, consistent with Maskobj
    '''
    lon = -7.0 # todo: get lon,lat from basins
    lat = 36.0
    jpk, _, _ = Maskobj.shape
    (nFrames, nSub, nCoast, jpk) = T.shape
    PRES = np.broadcast_to(Maskobj.zlevels,  (nFrames,nSub,nCoast,jpk))
    LON = np.broadcast_to(lon, (nFrames,nSub,nCoast,jpk))
    LAT = np.broadcast_to(lat, (nFrames,nSub,nCoast,jpk))
    SA = gsw.SA_from_SP(S, PRES,LON,LAT)
    CT = gsw.CT_from_pt(SA, T)
    RHO = gsw.rho(SA, CT, PRES)
    return RHO



if __name__=="__main__":
    maskfile="/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask_CMCC.nc"
    TheMask = Mask.from_file(maskfile, e3t_var_name="e3t_0")
    filename="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/AVE_PHYS/ave.20150116-12:00:00.phys.nc"
    rho =get_density(filename, TheMask)
    fileS = "/g100_scratch/userexternal/gbolzon0/V13C/GIB_profiles/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/vosaline.pkl"
    fileT = "/g100_scratch/userexternal/gbolzon0/V13C/GIB_profiles/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/votemper.pkl"
    TIMESERIES_T, TL = read_pickle_file(fileT)
    TIMESERIES_S, _ = read_pickle_file(fileS)
    rho = get_density_from_4D_profiles(TIMESERIES_T[:,:,:,:,0], TIMESERIES_S[:,:,:,:,0], TheMask)
    L = [rho, TL]
    import pickle
    with open("rho.pkl", "wb") as f:
        pickle.dump(L, f)
