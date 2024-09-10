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

class ncreader():
    def __init__(self, filename):
        ncIN = NC.netcdf_file(filename,"r")
        self.nTimes = ncIN.dimensions['time']
        self.nStats = ncIN.dimensions['nStats']
        self.nVar    = ncIN.dimensions['var']

        self.STATLIST=ncIN.statlist.split(",")
        self.VARLIST = ncIN.varlist.split(",")

        self.model = ncIN.variables['model'].data.copy()
        self.float   = ncIN.variables['float'   ].data.copy()
    def plotdata(self,var,stat):
        '''
        Arguments:
        * var  * string, modelvarname
        * stat * string  e.g. Int_0-200,Corr,DCM,z_01,Nit_1
        Returns:
        * Model * numpy array
        * Ref * numpy array
        '''
        
        ivar = self.VARLIST.index(var)
        istat = self.STATLIST.index(stat)
        Model = self.model[ivar,:,istat]
        Ref   = self.float[ivar,:,istat]
        return Model, Ref

if __name__=="__main__":
    A = ncreader("/pico/scratch/userexternal/lfeudale/validation/work/output/6901600.nc")
    model, biofloat =A.plotdata('P_l','Int_0-200')