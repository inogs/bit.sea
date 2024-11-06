from scipy.interpolate import griddata
import numpy as np
from bitsea.commons.utils import data_for_linear_interp
from scipy import interpolate

def shift(M2d,pos, verso):
    out = np.ones_like(M2d)*np.nan

    if (verso=='u'): out[:-pos,:] = M2d [pos:,:]
    if (verso=='d'): out[ pos:,:] = M2d[:-pos,:]
    if (verso=='l'): out[:,:-pos] = M2d [:,pos:]
    if (verso=='r'): out[:,pos:]  = M2d[:,:-pos]

    if (verso=="ul") : out[:-pos,:-pos] =  M2d [pos:,pos:]
    if (verso=="ur") : out[:-pos,pos: ] =  M2d [pos:,:-pos]
    if (verso=="dl") : out[ pos:,pos:]  =  M2d [:-pos,pos:]
    if (verso=="dr") : out[pos:,pos: ]  = M2d [:-pos,:-pos]
    return out



def SeaOverLand(M2d,n):
    jpj,jpi = M2d.shape
    a = np.zeros((8*n, jpj,jpi),np.float32)
    for i in range(n):
        a[i,:,:] = shift(M2d,i+1,'u')
    for i in range(n):
        a[i+n,:,:] = shift(M2d,i+1,'d')
    for i in range(n):
        a[i+2*n,:,:] = shift(M2d,i+1,'r')
    for i in range(n):
        a[i+3*n,:,:] = shift(M2d,i+1,'l')
    for i in range(n):
        a[i+4*n,:,:] = shift(M2d,i+1,'ur')
    for i in range(n):
        a[i+5*n,:,:] = shift(M2d,i+1,'dl')
    for i in range(n):
        a[i+6*n,:,:] = shift(M2d,i+1,'ul')
    for i in range(n):
        a[i+7*n,:,:] = shift(M2d,i+1,'dr')
    result = np.nanmean (a,axis=0)
    ii=~np.isnan(M2d)
    result[ii] = M2d[ii]
    return result

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
     
    points[:,0] = Mask1.xlevels[tmask1]
    points[:,1] = Mask1.ylevels[tmask1]
    values = Map2d[tmask1]
    MAP2d_nearest =griddata( points, values, (Mask2.xlevels, Mask2.ylevels), 'nearest', fill_value=np.nan)
    MAP2d_nearest[~tmask2]=np.nan
    return MAP2d_nearest


def interp_same_resolution(Mask1, Mask2, M3d):
    '''
    Performs nearest interpolation for masks very similar,
    differing only for sea/lands points.
    Arguments:
    * Mask1, Mask2 *   Mask objects
    * M3d          * 3d numpy array, consistent with Mask1

    Returns:
     * OUT * 3d numpy array, consistent with Mask2
    '''
    ii=Mask2.mask & (~Mask1.mask)
    K,J,I = np.nonzero(ii)
    OUT = M3d
    for k in np.unique(K):
        M2d=M3d[k,:,:]
        goods  = Mask1.mask[k,:,:]
        Jgoods, Igoods = np.nonzero(goods)
        nP = len(Jgoods)
        points = np.zeros((nP,2),dtype=np.float32)
        points[:,0] = Jgoods
        points[:,1] = Igoods
        values = M2d[goods]

        bool_mask2lands_on_k = ii[k,:,:]
        J,I = np.nonzero(bool_mask2lands_on_k)
        nP =bool_mask2lands_on_k.sum()
        xi = np.zeros((nP,2),dtype=np.float32)
        xi[:,0] = J
        xi[:,1] = I
        V = griddata(points, values, xi, "nearest")
        OUT[k,bool_mask2lands_on_k] = V
    return OUT

def space_interpolator_griddata(mask2, mask1, M3d):
    """Interpolates a 3d matrix using horizontal slices"""
    X1, Y1 = np.meshgrid(mask1.lon, mask1.lat)
    X2, Y2 = np.meshgrid(mask2.lon, mask2.lat)
    M3d_out = np.zeros((mask2.jpk, mask2.jpj, mask2.jpi), np.float32) * np.nan

    for k in range(mask2.jpk):
        jkb, jka, t_interp = data_for_linear_interp(
            mask1.zlevels, mask2.zlevels[k]
        )
        # indipendent from umask, vmask, or tmask
        tmask1 = (M3d[jkb, :, :] < 1.0e19) & ~np.isnan(M3d[jkb, :, :])
        # tmask1 = mask1.tmask[jkb,:,:]
        if tmask1.sum() == 0:
            print("All nans, return to upper layer ")
            jkb = jkb - 1
            # tmask1 = mask1.tmask[jkb,:,:]
            tmask1 = (M3d[jkb, :, :] < 1.0e19) & ~np.isnan(M3d[jkb, :, :])

        Map2d = M3d[jkb, :, :]
        nP = tmask1.sum()
        points = np.zeros((nP, 2), np.float32)
        points[:, 0] = X1[tmask1]
        points[:, 1] = Y1[tmask1]
        values = Map2d[tmask1]
        MAP2d_nearest = interpolate.griddata(
            points, values, (X2, Y2), "nearest", fill_value=np.nan
        )
        M3d_out[k, :, :] = MAP2d_nearest

    # M3d_out[~mask2.tmask] = np.nan # in order to avoid problems

    if np.isnan(M3d_out[mask2.mask]).any():
        print(
            "nans in space_interpolator_griddata:",
            np.isnan(M3d_out[mask2.mask]).sum(),
        )
        for k in range(mask2.jpk):
            a = M3d_out[k, :, :]
            lev_mask = mask2.mask[k, :, :]
            print(k, np.isnan(a[lev_mask]).sum())

    return M3d_out

def vertical_plane_interpolator(mask2, mask1, M2d, side):
    """Interpolates a 2d vertical matrix by using horizontal profiles
    M2d size is [jpk, LonSize, or LatSize]
    """

    if side in ["E", "W"]:
        M = np.zeros((mask2.jpk, mask2.jpj), dtype=np.float32)
        X1 = mask1.lat
        X2 = mask2.lat
    if side in ["N", "S"]:
        M = np.zeros((mask2.jpk, mask2.jpi), dtype=np.float32)
        X1 = mask1.lon
        X2 = mask2.lon

    if np.isnan(M2d).all():
        M[:, :] = np.nan
        return M

    for jk in range(mask2.jpk):
        jkb, jka, t_interp = data_for_linear_interp(
            mask1.zlevels, mask2.zlevels[jk]
        )
        Horizontal_Profile_b = M2d[jkb, :]
        Horizontal_Profile_a = M2d[jka, :]
        waterpoints_b = ~np.isnan(Horizontal_Profile_b)
        waterpoints_a = ~np.isnan(Horizontal_Profile_a)

        if waterpoints_b.any():
            Horiz_Profile_new_b = np.interp(
                X2, X1[waterpoints_b], Horizontal_Profile_b[waterpoints_b]
            )
        if waterpoints_a.any():
            Horiz_Profile_new_a = np.interp(
                X2, X1[waterpoints_a], Horizontal_Profile_b[waterpoints_a]
            )

        M[jk, :] = (
            Horiz_Profile_new_b * (1 - t_interp)
            + Horiz_Profile_new_a * t_interp
        )

        # Non lo calcola e si tiene il precedente, che ha sicuramente gia' calcolato
    return M


if __name__ == "__main__":
    from bitsea.commons import netcdf3    
    from bitsea.commons.mask import Mask
    from bitsea.commons.dataextractor import DataExtractor

    Mask24=Mask("/gpfs/scratch/userexternal/plazzari/eas_v12/eas_v12_6/wrkdir/MODEL/meshmask.nc",dzvarname="e3t_0")    
    Mask16=Mask('/gpfs/scratch/userexternal/gbolzon0/RA_COAST_03/wrkdir/MODEL/meshmask.nc')
    
    filename="/gpfs/scratch/userexternal/gbolzon0/RA_COAST_03/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/ave.20120816-12:00:00.P_l.nc"


    VAR = DataExtractor(Mask16,filename,"P_l").values
    VAR24 = surf_interp_2d(Mask16, Mask24, VAR[0,:,:])
    
    netcdf3.write_2d_file(VAR24, 'P_l', 'P_l_24.nc', Mask24, fillValue=1e+20)



