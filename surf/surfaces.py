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
                        break
                if not(np.isnan(matrixDCM[jj,ji])):
                    inddcm = profile>(matrixCM[jj,ji]-(matrixCM[jj,ji]-profile[0])/2.)
                    dzp = maskobj.e3t[:profile_len,jj,ji]
                    matrixDCMthick[jj,ji] = np.nansum(dzp[inddcm])
                    matrixIDCM[jj,ji] = np.nanargmin(np.abs(profile-matrixDCM[jj,ji]))

    matrixDCM[~tmask] = np.nan
    matrixCM[~tmask] = np.nan
    matrixIDCM[~tmask] = np.nan
    matrixDCMthick[~tmask] = np.nan
    return matrixDCM,matrixCM,matrixIDCM,matrixDCMthick

def MWB(chl,maskobj):
    '''
    Calculates Mixed Winter Bloom depth

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


def MWB2(chl,maskobj):
    '''
    Calculates Mixed Winter Bloom depth
    and the mean chl
    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    indlev = maskobj.getDepthIndex(300)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixMWB = np.zeros((jpj,jpi))
    matrixMWB[:,:] = np.nan
    matrixIMWB = np.zeros((jpj,jpi))
    matrixIMWB[:,:] = np.nan
    matrixI9MWB = np.zeros((jpj,jpi))
    matrixI9MWB[:,:] = np.nan
    matrixCMWB = np.zeros((jpj,jpi))
    matrixCMWB[:,:] = np.nan
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
                for iic, cc in enumerate(profile):
                    if (cc<=profile[0]*.9):
                        #matrixCMWB[jj,ji] = np.nanmean(profile[:iic])
                        matrixCMWB[jj,ji] = profile[0]
                        matrixI9MWB[jj,ji] = iic
                        break

    matrixMWB[~tmask] = np.nan
    matrixIMWB[~tmask] = np.nan
    matrixCMWB[~tmask] = np.nan
    matrixI9MWB[~tmask] = np.nan
    return matrixMWB,matrixIMWB,matrixCMWB,matrixI9MWB


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
    based on maximum gradient

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixNCL = np.zeros((jpj,jpi))
    matrixNCL[:,:] = np.nan
    matrixN = np.zeros((jpj,jpi))
    matrixN[:,:] = np.nan
    #matrixINCL = np.zeros((jpj,jpi))
    #matrixINCL[:,:] = np.nan
    matrixNCLslope = np.zeros((jpj,jpi))
    matrixNCLslope[:,:] = np.nan
    matrixNCLthick = np.zeros((jpj,jpi))
    matrixNCLthick[:,:] = np.nan
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
                pslope = dN[iid]
                iip = dN>(.75*pslope)
                iia = np.argwhere(np.diff(iip))[:,0]
                if (len(iia))>1:
                    iit = np.zeros(len(dN),dtype=bool)
                    for ind in range(len(iia)-1):
                        if iid in range(iia[ind]+1,iia[ind+1]+1):
                            iit[iia[ind]+1:iia[ind+1]+1] = True
                            break
                    dzp = maskobj.e3t[:profile_len-1,jj,ji]
                    matrixNCLthick[jj,ji] = np.nansum(dzp[iit])
                    matrixNCLslope[jj,ji] = np.nanmean(dN[iit])

    
    matrixNCL[~tmask] = np.nan
    matrixN[~tmask] = np.nan
    #matrixINCL[~tmask] = np.nan
    matrixNCLslope[~tmask] = np.nan
    #return matrixNCL,matrixN,matrixINCL,matrixNCLslope
    return matrixNCL,matrixN,matrixNCLslope,matrixNCLthick



def NLAYERS(nut,maskobj):
    '''
    Calculates layers above and below nutricline and their average concentration

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixSURF = np.zeros((jpj,jpi))
    matrixSURF[:,:] = np.nan
    matrixNSURF = np.zeros((jpj,jpi))
    matrixNSURF[:,:] = np.nan
    matrixBOTT = np.zeros((jpj,jpi))
    matrixBOTT[:,:] = np.nan
    matrixNBOTT = np.zeros((jpj,jpi))
    matrixNBOTT[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = bathy_len
                profile = nut[:profile_len,jj,ji]
                depths = maskobj.zlevels[:profile_len]
                dN = np.diff(profile)/np.diff(depths)
                iid = np.nanargmax(dN)
                pslope = dN[iid]
                dNsurf = np.diff(profile[:iid])/np.diff(depths[:iid])
                iis = np.abs(dNsurf)<.1*pslope
                if np.sum(iis)>1:
                    matrixSURF[jj,ji] = depths[:iid-1][iis][-1]
                    e3tp = maskobj.e3t[:iid-1,jj,ji][iis]
                    profpp = profile[:iid-1][iis]
                    matrixNSURF[jj,ji] = np.nansum(profpp*e3tp)/np.nansum(e3tp)
                dNbott = np.diff(profile[iid:])/np.diff(depths[iid:])
                iib = np.abs(dNbott)<.1*pslope
                if np.sum(iib)>1:
                    matrixBOTT[jj,ji] = depths[iid:-1][iib][0]
                    e3tp = maskobj.e3t[iid:profile_len-1,jj,ji][iib]
                    profpp = profile[iid:-1][iib]
                    matrixNBOTT[jj,ji] = np.nansum(profpp*e3tp)/np.nansum(e3tp)

    matrixSURF[~tmask] = np.nan
    matrixBOTT[~tmask] = np.nan
    matrixNSURF[~tmask] = np.nan
    matrixNBOTT[~tmask] = np.nan
    return matrixSURF,matrixBOTT,matrixNSURF,matrixNBOTT
    #return matrixSURF,matrixBOTT

def PICNCL_dz_max(den,maskobj):
    '''
    Calculates picnocline depth

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixPCL = np.zeros((jpj,jpi))
    matrixPCL[:,:] = np.nan
    matrixD = np.zeros((jpj,jpi))
    matrixD[:,:] = np.nan
    #matrixIPCL = np.zeros((jpj,jpi))
    #matrixIPCL[:,:] = np.nan
    matrixPCLslope = np.zeros((jpj,jpi))
    matrixPCLslope[:,:] = np.nan
    matrixPCLthick = np.zeros((jpj,jpi))
    matrixPCLthick[:,:] = np.nan
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = bathy_len
                profile = den[:profile_len,jj,ji]
                depths = maskobj.zlevels[:profile_len]
                dN = np.diff(profile)/np.diff(depths)
                dN[depths[:-1]<=5] = np.nan
                iid = np.nanargmax(dN)
                PCL = depths[iid]
                DEN = profile[iid]
                matrixPCL[jj,ji] = PCL
                matrixD[jj,ji] = DEN
                #matrixIPCL[jj,ji] = iid
                pslope = dN[iid]
                iip = dN>(.75*pslope)
                iia = np.argwhere(np.diff(iip))[:,0]
                if (len(iia))>1:
                    iit = np.zeros(len(dN),dtype=bool)
                    for ind in range(len(iia)-1):
                        if iid in range(iia[ind]+1,iia[ind+1]+1):
                            iit[iia[ind]+1:iia[ind+1]+1] = True
                            break
                    dzp = maskobj.e3t[:profile_len-1,0,0]
                    matrixPCLthick[jj,ji] = np.nansum(dzp[iit])
                    matrixPCLslope[jj,ji] = np.nanmean(dN[iit])


    matrixPCL[~tmask] = np.nan
    matrixD[~tmask] = np.nan
    #matrixIPCL[~tmask] = np.nan
    matrixPCLslope[~tmask] = np.nan
    #return matrixPCL,matrixD,matrixIPCL,matrixPCLslope
    return matrixPCL,matrixD,matrixPCLslope,matrixPCLthick

def DIFFST(dff,maskobj):
    '''
    Calculates some statistics on diffusion 

    '''
    _,jpj,jpi = maskobj.shape
    tmask = maskobj.mask_at_level(200)
    # indlev = maskobj.getDepthIndex(0)
    DEPTHS = maskobj.bathymetry_in_cells()
    matrixMAX = np.zeros((jpj,jpi))
    matrixMAX[:,:] = np.nan
    matrixMEAN = np.zeros((jpj,jpi))
    matrixMEAN[:,:] = np.nan
    matrixDLOW = np.zeros((jpj,jpi))
    matrixDLOW[:,:] = np.nan
    i200 = maskobj.getDepthIndex(200)
    for jj in range(jpj):
        for ji in range(jpi):
            if tmask[jj,ji]:
                bathy_len = DEPTHS[jj,ji]
                profile_len = bathy_len
                profile = dff[:profile_len,jj,ji]
                profile[profile<1.e-3] = np.nan
                depths = maskobj.zlevels[:profile_len]
                idd = np.min([i200,profile_len])
                if not np.all(np.isnan(profile[:idd])):
                    imax = np.nanargmax(profile[:idd])
                    i10 = np.where(np.isnan(profile[:idd])==False)[0][-1]
                    matrixMAX[jj,ji] = profile[imax]
                    matrixDLOW[jj,ji] = depths[i10]
                    matrixMEAN[jj,ji] = np.nanmean(profile[:idd])

    matrixMAX[~tmask] = np.nan
    matrixMEAN[~tmask] = np.nan
    matrixDLOW[~tmask] = np.nan
    return matrixMAX,matrixMEAN,matrixDLOW


