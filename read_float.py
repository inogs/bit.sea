import numpy as np
import scipy.io.netcdf as NC
import pylab as pl
import sys
import datetime


VARLIST=['DOXY','NITRATE','CHLA', 'PRES','PSAL','TEMP']



    
    
def searchVariable_on_parameters(filename,var):
    '''
    returns index of profile on which variable has to be read, 
    -1 if fails (PARAMETERS variable is blank)
    
    '''
    
    ncIN = NC.netcdf_file(filename,'r')
    N_PROF   =ncIN.dimensions['N_PROF']
    N_PARAM  =ncIN.dimensions['N_PARAM']
    PARAMETER=ncIN.variables['PARAMETER'].data.copy()
    ncIN.close()
    if PARAMETER.tostring().rstrip() == '':
        iProf      = -1
        return iProf
    else:
        for iprof in range(N_PROF):
            for iparam in range(N_PARAM):
                s=PARAMETER[iprof,0,iparam,:].tostring().rstrip()
                if s==var:
                    return iprof        



def fillnan(ncObj,var):
    varObj = ncObj.variables[var]
    fillvalue  =varObj._FillValue
    M     = varObj.data.copy()
    M[M==fillvalue] = np.NaN;
    return M


def merge_profile_with_adjusted(profile, profile_adj):
    N_LEV = profile.size
    res = np.zeros_like(profile)
    if np.isnan(profile_adj).all():
        res = profile
    else:
        for k in range(N_LEV):
            if np.isnan(profile_adj[k]):
                res[k] = profile[k]
            else:
                res[k] = profile_adj[k]
    return res    

def merge_var_with_adjusted(ncObj,var):
    N_PROF= ncObj.dimensions['N_PROF']     
    M     = fillnan(ncObj, var)
    M_ADJ = fillnan(ncObj, var + "_ADJUSTED")
    M_RES = M
    for iprof in range(N_PROF):
        M_RES[iprof,:] = merge_profile_with_adjusted(M[iprof,:], M_ADJ[iprof,:]) 

    return M_RES
  
    

INDEX_FILE=np.loadtxt("index.txt",dtype=mydtype, delimiter=",")
nFiles=INDEX_FILE.size
var = 'NITRATE' 
   
for iFile in range(nFiles):
    filename = INDEX_FILE['file_name'][iFile]
    available_params = INDEX_FILE['parameters'][iFile]     
    if var in available_params:     
        iProf = searchVariable_on_parameters(filename, var)
        print iProf
        ncIN=NC.netcdf_file(filename,'r')
        
        if iProf== -1 :
            M = merge_var_with_adjusted(ncIN, var)
            N_PROF= ncIN.dimensions['N_PROF']
            N_MEAS = np.zeros((N_PROF),np.int32)
            for iprof in range(N_PROF):
                N_MEAS[iprof] = (~np.isnan(M[iprof,:])).sum()
            #per il momento cerco il massimo
            iProf = N_MEAS.argmax()
            print "iprof new value", iProf
            
        M = merge_var_with_adjusted(ncIN, var)
        PRES    = merge_var_with_adjusted(ncIN, 'PRES')
        Profile =  M[iProf,:]
        Pres    =  PRES[iProf,:]                            
        ncIN.close()
        
        #pl.figure(); pl.plot(Profile,Pres); pl.gca().invert_yaxis(); pl.show(block=False)
            
            