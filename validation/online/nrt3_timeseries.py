from commons.timeseries import TimeSeries
import glob
import numpy as np
import scipy.io.netcdf as NC
import datetime

class timelistcontainer():
    def __init__(self,Ti,ARCHIVE_DIR,postfix_dir=""):
        
        self.timelist=[]
        self.filelist=[]        
        self.bias = None
        self.number=None
        self.rmse  =None

        TS=TimeSeries(Ti,archive_dir=ARCHIVE_DIR, postfix_dir=postfix_dir, glob_pattern="BioFloat_Weekly_validation")
        search_paths = TS.get_runs([2])
        for _, directory in search_paths:
            dirlist=glob.glob(directory + "BioFloat_Weekly_validation*")
            if len(dirlist)> 0: 
                self.filelist.append(dirlist[0])
        self.nFrames = len(self.filelist)
        self.read_basic_info(self.filelist[0])
        self.readfiles()

    def read_basic_info(self,filename):
        ncIN = NC.netcdf_file(filename,'r')
        self.nVAR = ncIN.dimensions['var']
        self.nSUB = ncIN.dimensions['sub']
        self.nDEPTH  = ncIN.dimensions['depth']
        self.SUBLIST = ncIN.sublist.split(",")
        self.LAYERLIST=ncIN.layerlist.split(",")
        self.VARLIST = ncIN.varlist.split(",")
        
        ncIN.close()        

    def read_validation_file(self,filename):
        ncIN =NC.netcdf_file(filename,'r')
        NPOINTS=ncIN.variables['npoints'].data.copy()
        BIAS=ncIN.variables['bias'].data.copy()
        RMSE=ncIN.variables['rmse'].data.copy()
        return NPOINTS, BIAS, RMSE

    def readfiles(self):
        
        BIAS    =np.zeros((self.nFrames,self.nVAR, self.nSUB,self.nDEPTH))
        RMSE    =np.zeros((self.nFrames,self.nVAR, self.nSUB,self.nDEPTH))
        NUMBER  =np.zeros((self.nFrames,self.nVAR, self.nSUB,self.nDEPTH))

        
        for iFrame, filename in enumerate(self.filelist):
            time = datetime.datetime.strptime(filename[-11:],'%Y%m%d.nc')
            self.timelist.append(time)
            number, bias, rmse = self.read_validation_file(filename)
            NUMBER[iFrame,:,:] = number
            BIAS[iFrame,:,:]   = bias
            RMSE[iFrame,:,:]   = rmse

        self.bias = BIAS
        self.number=NUMBER
        self.rmse = RMSE

    def plotdata(self,VAR, var,sub,depth):
        '''
        VAR must be a 3D field of this class, such as bias, number, ...
        '''
        ivar = self.VARLIST.index(var)
        isub = self.SUBLIST.index(sub)
        idepth= self.LAYERLIST.index(depth)
        return self.timelist, VAR[:,ivar, isub,idepth]



if __name__ == '__main__':
    from commons.time_interval import TimeInterval
    TI = TimeInterval("20150608","20160101","%Y%m%d")
    ARCHIVEDIR="/pico/scratch/userexternal/gbolzon0/NRT/V4/NRT3_outputs/"
    A = timelistcontainer(TI,ARCHIVEDIR)
    A.plotdata(A.bias, 'P_i','tyr', '10-30m')
