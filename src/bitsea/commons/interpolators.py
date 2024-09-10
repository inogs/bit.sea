from scipy.interpolate import griddata
import numpy as np
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



