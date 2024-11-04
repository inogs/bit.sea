import numpy as np
from bitsea.Sat import SatManager
from scipy import interpolate

def get_2_indices_for_slicing(array,MinValue,MaxValue, istart):
    n = len(array)
    for i in range(istart, n):
        if array[i]>= MinValue:
            MinIndex=i
            break
    for i in range(MinIndex,n):
        if array[i]>= MaxValue:
            MaxIndex=i
            break
        MaxIndex = i

    return MinIndex, MaxIndex

def array_of_indices_for_slicing(xcoarse, xfine):
    '''
    Returns two arrays, having the same size of xcoarse
    I_START, I_END
    
    '''
    
    jpi = len(xcoarse)
    deltax = np.diff(xcoarse).mean() # 1./16
    I_START = np.zeros((jpi,),int)
    I_END   = np.zeros((jpi,),int)
    istart = 0
    for ji in range(jpi):
        centrocella = xcoarse[ji]
        bordoW = centrocella-deltax/2
        bordoE = centrocella+deltax/2
        try:
            istart, iend = get_2_indices_for_slicing(xfine, bordoW, bordoE, istart)
        except:
            print("*******************************************************")
            print("Warning: out of bounds.")
            print("We will use the previous istart, iend values")
            print("COARSE GRID: Cell center, %f, West bound, %f East bound %f "  %(centrocella, bordoW, bordoE))
            print("FINE  GRID:  ends at: ", xfine.max())
            print("*******************************************************")
        I_START[ji] = istart
        I_END[ji]   = iend
    return I_START, I_END


def interp_2d_by_cells_slices(Mfine, Maskout, I_START, I_END, J_START, J_END, fillValue=-999.0, min_cov=0.0, ave_func=SatManager.mean):
    '''
    Interpolates data from a fine mesh to a coarser one.
    
    Arguments:
    * Mfine   * the 2d matrix to interpolate
    * Maskout * mask object of Coarse mask, to get tmask
    * I_START, I_END, J_START, J_END * array of integers defining box
    * min_cov  * float, between 0 and 1, minimum of accepted coverage.
    * ave_func * a function to get an average

    Returns:
    * OUT *  2d matrix interpolated on Maskout
    *  NP *  Number of used points, 2d matrix of integers
    '''
    _, jpj, jpi = Maskout.shape
    tmask = Maskout.mask_at_level(0)
    OUT = np.ones((jpj,jpi),np.float32)*fillValue
    NP = np.ones((jpj,jpi),np.int32)
    
    for ji in range(jpi):
        istart = I_START[ji]
        i_end  =   I_END[ji]
        # istart == iend if Coarse grid is out of Fine grid
        if istart == i_end : i_end=istart+1 #continue
        for jj in range(jpj):
            jstart = J_START[jj]
            j_end  =   J_END[jj]
            if j_end==jstart: j_end=jstart+1
            if tmask[jj,ji]:
                localcell = Mfine[jstart:j_end, istart:i_end]
                goods = localcell>0
                nPoints = goods.sum()
                NP[jj,ji] = nPoints
                if np.any(goods):
                    coverage = float(nPoints)/localcell.size
                    if (coverage >= min_cov):
                        OUT[jj,ji] = ave_func(localcell[goods])
    return OUT, NP

def nearest(Min, MaskIn, Maskout):
    M = Min.copy()
    x = MaskIn.lon
    y = MaskIn.lat
    goods = M > 0
    M[ ~ goods] = np.nan
    X,Y = np.meshgrid(x,y)
    nPoints = X.size
    Xpoints = np.zeros((nPoints,2),float)
    Xpoints[:,0] = X.ravel()
    Xpoints[:,1] = Y.ravel()
    interp = interpolate.NearestNDInterpolator(Xpoints, M.ravel())
    _,jpj,jpi = Maskout.shape
    Z = interp(Maskout.xlevels, Maskout.ylevels)


    NP = np.ones((jpj,jpi),int)
    NP[np.isnan(Z)] = 0
    return Z, NP

if __name__== "__main__":

    from bitsea.commons.mask import Mask
    from bitsea.Sat import SatManager as Sat
    TheMask=Mask('/leonardo_scratch/large/userexternal/gbolzon0/MIT/V1/devel/wrkdir/BC_IC/mask.nc')

    jpk,jpj,jpi = TheMask.shape
    x = TheMask.lon
    y = TheMask.lat


    x1km = Sat.masks.KD490mesh.lon
    y1km = Sat.masks.KD490mesh.lat

    I_START, I_END = array_of_indices_for_slicing(x, x1km)
    J_START, J_END = array_of_indices_for_slicing(y, y1km)


    inputfile="/leonardo_scratch/large/userexternal/gbolzon0/ONLINE//SAT/CHL/DT/WEEKLY_1_1km/20241021_cmems_obs-oc_med_bgc-plankton_myint_l3-multi-1km_P1D.nc"
    CHL = Sat.readfromfile(inputfile,'CHL')

    #Mout, Usedpoints =interp_2d_by_cells_slices(CHL, TheMask, I_START, I_END, J_START, J_END)
    Mout, Usedpoints = nearest(CHL, Sat.masks.KD490mesh, TheMask )
    Sat.dumpGenericfile('test.nc', Mout, 'CHL', Sat.masks.CadeauMesh)
        

            
    
