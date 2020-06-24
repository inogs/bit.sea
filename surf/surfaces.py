import numpy as np
from commons.mask import Mask
import sys
def mld(temperature,maskobj):
    ''' Calculates of Mixed Layer Depth based on temperature 
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
    Rough calculates of Deep Chlorophyll Maximum '''

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
    Improved calculates of Deep Chlorophyll Maximum
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
    Calculates of Deep Chlorophyll Maximum
    Uses 1-st and 2-nd derivative to find the maximum
    Returns the cchlorophyll concentration at DCM also

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    indlev = maskobj.getDepthIndex(200)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixDCM = np.zeros((jpj,jpi))
    matrixDCM[:,:] = np.nan
    matrixCM = np.zeros((jpj,jpi))
    matrixCM[:,:] = np.nan
    matrixIDCM = np.zeros((jpj,jpi))
    matrixIDCM[:,:] = np.nan
    matrixDCMthick = np.zeros((jpj,jpi))
    matrixDCMthick[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = min(bathy_len,indlev+1)
                maskprof = chl[:profile_len,jj,ji]>0.1
                profile = chl[:profile_len,jj,ji]
                depths = maskobj.zlevels[:profile_len]
                profile_filt = profile[maskprof]
                depths_filt = depths[maskprof]
                profile_rev = profile_filt[::-1]# from bottom
                depths_rev = depths_filt[::-1]
                d1 = np.diff(profile_rev,1)
                d2 = np.diff(profile_rev,2)
                for iid, dd in enumerate(d1):
                    if (iid>0) and (dd<0) and (d2[iid-1]<0):
                        max_cand = profile_rev[iid-1:iid+2]
                        d_max_cand = depths_rev[iid-1:iid+2]
                        indmax = np.argmax(max_cand)
                        CM = max_cand[indmax]
                        DCM = d_max_cand[indmax]
                        matrixDCM[jj,ji] = DCM
                        matrixCM[jj,ji] = CM
                        matrixIDCM[jj,ji] = profile_len-(iid+indmax)
                        break
                if not(np.isnan(matrixDCM[jj,ji])):
                    inddcm = profile>(matrixCM[jj,ji]-(matrixCM[jj,ji]-profile[0])/2.)
                    dzp = maskobj.e3t[:profile_len,jj,ji]
                    matrixDCMthick[jj,ji] = np.nansum(dzp[inddcm])

    matrixDCM[~tmask] = np.nan
    matrixCM[~tmask] = np.nan
    matrixIDCM[~tmask] = np.nan
    matrixDCMthick[~tmask] = np.nan
    return matrixDCM,matrixCM,matrixIDCM,matrixDCMthick

def MWB(chl,maskobj):
    '''
    Calculates of Mixed Winter Bloom depth
    Uses 1-st and 2-nd derivative to find the maximum
    Returns the chlorophyll concentration at DCM also

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    indlev = maskobj.getDepthIndex(300)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixMWB = np.zeros((jpj,jpi))
    matrixMWB[:,:] = np.nan
    matrixIMWB = np.zeros((jpj,jpi))
    matrixIMWB[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = min(bathy_len,indlev+1)
                profile = chl[:profile_len,jj,ji]
                depths = maskobj.zlevels[:profile_len]
                for iid, dd in enumerate(profile):
                    if (dd<=profile[0]*.1):
                        MWB = depths[iid]
                        matrixMWB[jj,ji] = MWB
                        matrixIMWB[jj,ji] = iid
                        break
                
    
    matrixMWB[~tmask] = np.nan
    matrixIMWB[~tmask] = np.nan
    return matrixMWB,matrixIMWB


def NITRCL(nit,maskobj):
    '''
    Calculates of nitracline depth
    based on criteria p>=2

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    # indlev = maskobj.getDepthIndex(0)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixNCL = np.zeros((jpj,jpi))
    matrixNCL[:,:] = np.nan
    matrixINCL = np.zeros((jpj,jpi))
    matrixINCL[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = bathy_len
                profile = nit[:profile_len,jj,ji]
                depths = maskobj.zlevels[:profile_len]
                for iid, dd in enumerate(profile):
                    if (dd>=2) and (depths[iid]>30):
                        NCL = depths[iid]
                        matrixNCL[jj,ji] = NCL
                        matrixINCL[jj,ji] = iid
                        break
                
    
    matrixNCL[~tmask] = np.nan
    matrixINCL[~tmask] = np.nan
    return matrixNCL,matrixINCL


def NUTRCL_dz_max(nut,maskobj):
    '''
    Calculates of nitracline depth
    based on criteria p>=2

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixNCL = np.zeros((jpj,jpi))
    matrixNCL[:,:] = np.nan
    matrixN = np.zeros((jpj,jpi))
    matrixN[:,:] = np.nan
    matrixINCL = np.zeros((jpj,jpi))
    matrixINCL[:,:] = np.nan
    matrixNCLslope = np.zeros((jpj,jpi))
    matrixNCLslope[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = bathy_len
                profile = nut[:profile_len,jj,ji]
                depths = maskobj.zlevels[:profile_len]
                dN = np.diff(profile)/np.diff(depths)
                dN[depths[:-1]<=30] = np.nan
                iid = np.nanargmax(dN)
                NCL = depths[iid]
                NUT = profile[iid]
                matrixNCL[jj,ji] = NCL
                matrixN[jj,ji] = NUT
                matrixINCL[jj,ji] = iid
                matrixNCLslope[jj,ji] = np.nanmean(dN[iid-2:iid+2])
 
    
    matrixNCL[~tmask] = np.nan
    matrixN[~tmask] = np.nan
    matrixINCL[~tmask] = np.nan
    matrixNCLslope[~tmask] = np.nan
    return matrixNCL,matrixN,matrixINCL,matrixNCLslope

