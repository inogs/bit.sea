import numpy as np
from commons.mask import Mask

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
    Rough calculation of Deep Clorophyll Maximum '''

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
                end_index=mindiff2+3
                if end_index>profile_len:
                    end_index=profile_len
                
                local_ind = np.argmax(profile[mindiff2:end_index])
                matrix_2D[jj,ji] = profile_len - (local_ind + mindiff2)
                
    out = maskobj.zlevels[matrix_2D]
    out[~tmask] = 1.e+20

