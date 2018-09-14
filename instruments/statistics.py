import numpy as np
def mean_profile(ProfileList, var, depth):
    '''
    Calculates the mean profile of Profilelist over a provided depth
    Arguments:
    * ProfileList * a list of Profile Objects
    * var         * string, the var name
    * depth       * numpy array, the depth of the output

    Returns:
    * MEAN * numpy array corresponding to depth where profiles of measurements are interpolated
    '''
    nProfiles = len(ProfileList) 
    nLev = len(depth)
    
    M = np.zeros((nProfiles, nLev), np.float32)*np.nan
    MEAN = np.zeros((nLev,),np.float32) * np.nan
    
    for ip, p in enumerate(ProfileList):
        Pres, Profile, Qc = p.read(var)
        M[ip,:] = np.interp(depth, Pres, Profile, right=np.nan)
    for iLev in range(nLev):
        level= M[:,iLev]
        good = ~np.isnan(level)
        if good.sum() > 0:
            MEAN[iLev] = level[good].mean()
    return MEAN

