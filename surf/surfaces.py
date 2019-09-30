import numpy as np
from commons.mask import Mask
import sys
def mld(temperature,maskobj):
    ''' Calculation of Mixed Layer Depth based on temperature 
    mld is defined as 
     '''
    jpk,jpj,jpi=maskobj.shape
    tmask=maskobj.mask_at_level(0)
    DEPTHS=maskobj.bathymetry_in_cells()
    matrix_2D = np.ones((jpj,jpi),np.int)*(jpk-1)
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                depth_cell=DEPTHS[jj,ji]
                absdiff_array =  abs(temperature[:depth_cell,jj,ji]- temperature[3,jj,ji])
                for k,absdiff in enumerate(absdiff_array):
                    if absdiff > 0.1:
                        break
                matrix_2D[jj,ji] = k

    Mixed_Layer_Depth=maskobj._zlevels[matrix_2D]
    max_value=500
    ii = Mixed_Layer_Depth > max_value
    Mixed_Layer_Depth[ii] = max_value
    Mixed_Layer_Depth[~tmask] = 1.e+20
    return Mixed_Layer_Depth

def dcm(chl,maskobj):
    '''
    Rough calculation of Deep Chlorophyll Maximum '''

    jpk,jpj,jpi=maskobj.shape
    tmask   = maskobj.mask
    tmask2d = tmask[0,:,:]
    chl[~tmask]=0.0
    #indexes==chl.argmax(axis=0)
    chl=chl[-1::-1,:,:]
    reverted_indexes=chl.argmax(axis=0)
    indexes=jpk-1-reverted_indexes
    out = maskobj.zlevels[indexes]
    out[~tmask2d] = 1.e+20
    return out

def DCM(chl,maskobj):
    '''
    Improved calculation of Deep Chlorophyll Maximum
    Uses 2-nd derivative to find the first maximum, starting from bottom.

    '''
    jpk,jpj,jpi=maskobj.shape
    tmask=maskobj.mask_at_level(0)
    DEPTHS=maskobj.bathymetry_in_cells()
    matrix_2D = np.ones((jpj,jpi),np.int)*(jpk-1)
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                profile_len = DEPTHS[jj,ji]
                profile = chl[:profile_len,jj,ji]
                profile = profile[-1::-1]# from bottom
                d2 = np.diff(profile,2)
                mindiff2=np.argmin(d2)
                if d2[mindiff2] >  0: # non trovo un massimo
                     matrix_2d = profile_len - np.argmax(profile)
                else:
                    end_index=mindiff2+3
                    if end_index>profile_len:
                        end_index=profile_len
                
                    local_ind = np.argmax(profile[mindiff2:end_index])
                    matrix_2D[jj,ji] = profile_len - (local_ind + mindiff2)
                
    out = maskobj.zlevels[matrix_2D]
    out[~tmask] = 1.e+20
    return out

def DCM2(chl,maskobj):
    '''
    Calculation of Deep Chlorophyll Maximum
    Uses 1-st and 2-nd derivative to find the maximum
    Returns the cchlorophyll concentration at DCM also

    '''
    _,jpj,jpi=maskobj.shape
    tmask=maskobj.mask_at_level(0)
    DEPTHS=maskobj.bathymetry_in_cells()
    matrixDCM = np.ones((jpj,jpi),np.int)
    matrixDCM[:,:] = np.nan
    matrixCM = np.ones((jpj,jpi),np.int)
    matrixCM[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                CM = np.nan
                DCM = np.nan
                profile_len = DEPTHS[jj,ji]
                maskprof = chl[:profile_len,jj,ji]>0.1
                profile = chl[:profile_len,jj,ji]
                depths = moskobj.zlevels[:profile]
                profile_filt = profile[maskprof]
                depths_filt = depths[maskprof]
                profile_rev = profile_filt[::-1]# from bottom
                depths_rev = depths_filt[::-1]
                d1 = np.diff(profile_rev,1)
                d2 = np.diff(profile_rev,2)
                # maskd1sign = np.sign(d1)>=0
                # mindiff2=np.argmin(d2)
                for iid, dd in range(d1):
                    if (iid>0) and (dd<0) and(d2[iid-1]<0):
                        max_cand = profile_rev[iid-1:iid+2]
                        d_max_cand = depths_rev[iid-1:iid+2]
                        indmax = np.argmax(max_cand)
                        CM = max_cand[indmax]
                        DCM = d_max_cand[indmax]
                        matrixDCM[jj,ji] = DCM
                        matrixCM[jj,ji] = CM
                
    
    matrixDCM = [~tmask] = 1.e+20
    matrixCM = [~tmask] = 1.e+20
    return matrixDCM,matrixCM

