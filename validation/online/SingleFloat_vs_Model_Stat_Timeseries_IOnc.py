import scipy.io.netcdf as NC
import numpy as np
def dumpfile(filename,M_FLOAT, M_MODEL, STATLIST):
    '''
    Arguments:
    * filename * string
    * M_FLOAT  * [nVar, nTimes,nStats] numpy array
    * M_MODEL  * [nVar, nTimes,nStats] numpy array
    * STATLIST * list of strings
    '''
    
    nFrames, nStats = M_FLOAT.shape
    ncOUT = NC.netcdf_file(filename,'w') #manca l'array times
    ncOUT.createDimension('time', nFrames)
    ncOUT.createDimension('nStats', nStats)
    
    s='';
    for stat in STATLIST: s =s+stat + ","
    setattr(ncOUT,'statlist',s[:-1])
    
    
    ncvar=ncOUT.createVariable('float', 'f', ('time', 'nStats'))
    ncvar[:] = M_FLOAT
    ncvar=ncOUT.createVariable('model', 'f', ('time','nStats'))
    ncvar[:] = M_MODEL
    
    ncOUT.close()

class ncreader():
    def __init__(self, filename):
        ncIN = NC.netcdf_file(filename,"r")
        self.nTimes = ncIN.dimensions['time']
        self.nStats = ncIN.dimensions['nStats']

        self.STATLIST=ncIN.statlist.split(",")

        self.model = ncIN.variables['model'].data.copy()
        self.float   = ncIN.variables['float'   ].data.copy()
        ncIN.close()
    def plotdata(self,var,stat, only_good=False):
        '''
        Arguments:
        * var  * string, modelvarname
        * stat * string  e.g. Int_0-200,Corr,DCM,z_01,Nit_1
        * only_good * boolean if True, returns data without nans, else the whole array 
        Returns:
        * Model * numpy array
        * Ref * numpy array
        '''
        
        istat = self.STATLIST.index(stat)
        Model = self.model[:,istat]
        Ref   = self.float[:,istat]
        if only_good==False:
            return Model, Ref
        else:
            ii= ~np.isnan(self.float[:,0])
            return Model[ii], Ref[ii]

if __name__=="__main__":
    A = ncreader("/pico/scratch/userexternal/lfeudale/validation/work/output/6901600.nc")
    model, biofloat =A.plotdata('P_l','Int_0-200')