import scipy.io.netcdf as NC
def dumpfile(filename,M_FLOAT, M_MODEL, VARLIST, STATLIST):
    '''
    Arguments:
    * filename * string
    * M_FLOAT  * [nVar, nTimes,nStats] numpy array
    * M_MODEL  * [nVar, nTimes,nStats] numpy array
    * VARLIST  * list of strings
    * STATLIST * list of strings
    '''
    
    nVar, nFrames, nStats = M_FLOAT.shape
    ncOUT = NC.netcdf_file(filename,'w') #manca l'array times
    ncOUT.createDimension('time', nFrames)
    ncOUT.createDimension('var', nVar)
    ncOUT.createDimension('nStats', nStats)
    
    s=''
    for var in VARLIST: s= s+var + ","
    setattr(ncOUT, 'varlist',s[:-1])
    s='';
    for stat in STATLIST: s =s+stat + ","
    setattr(ncOUT,'statlist',s[:-1])
    
    
    ncvar=ncOUT.createVariable('float', 'f', ('var','time', 'nStats'))
    ncvar[:] = M_FLOAT
    ncvar=ncOUT.createVariable('model', 'f', ('var','time','nStats'))
    ncvar[:] = M_MODEL
    
    ncOUT.close()
