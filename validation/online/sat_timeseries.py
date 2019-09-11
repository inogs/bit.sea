from commons.timeseries import TimeSeries
from commons.Timelist import TimeList
import glob
import numpy as np
import scipy.io.netcdf as NC
import datetime


class timelistcontainer():
    def __init__(self,Ti,ARCHIVE_DIR,searchstring,prefix="Validation_f1_20190507_on_daily_Sat"):
        '''
        timelistcontainer is a reader class for files
        produced by SatValidation, in operational chain or similar sequential runs
        Arguments:
        * Ti * a TimeInterval
        * ARCHIVE_DIR * string indicating path of the chain archive directory e.g. /pico/home/usera07ogs/a07ogs00/OPA/V2C/archive/
        * search_type * string, having one of these values: f0, f1, f2, a0, a1, a2
        * postfix_dir * string to pass to TimeSeries objects
        '''


        self.timelist=[]
        self.filelist=[]
        self.bias = None
        self.number=None
        self.rmse  =None
        self.model =None
        self.sat  = None
        self.rmselog = None
        self.biaslog = None
        self.search_type=searchstring
        TL = TimeList.fromfilenames(Ti, ARCHIVE_DIR, searchstring, prefix=prefix, dateformat="%Y%m%d")

        self.filelist = TL.filelist
        self.nFrames = len(self.filelist)
        self.read_basic_info(self.filelist[0])
        self.readfiles()

    def initold(self,Ti,ARCHIVE_DIR,search_type,postfix_dir=""):
        '''
        timelistcontainer is a reader class for files
        produced by SatValidation, in operational chain or similar sequential runs
        Arguments:
        * Ti * a TimeInterval 
        * ARCHIVE_DIR * string indicating path of the chain archive directory e.g. /pico/home/usera07ogs/a07ogs00/OPA/V2C/archive/
        * search_type * string, having one of these values: f0, f1, f2, a0, a1, a2
        * postfix_dir * string to pass to TimeSeries objects
        '''


        self.timelist=[]
        self.filelist=[]        
        self.bias = None
        self.number=None
        self.rmse  =None
        self.model =None
        self.sat  = None
        self.rmselog = None
        self.biaslog = None
        self.search_type=search_type
        TS=TimeSeries(Ti,archive_dir=ARCHIVE_DIR, postfix_dir=postfix_dir, glob_pattern="Validation_f0")
        search_paths = TS.get_runs([2])
        for _, directory in search_paths:
            dirlist=glob.glob(directory + "Validation_" + search_type + "*")
            if len(dirlist)> 0: 
                self.filelist.append(dirlist[0])
        self.nFrames = len(self.filelist)
        self.read_basic_info(self.filelist[0])
        self.readfiles()
    
    def read_basic_info(self,filename):
        ncIN = NC.netcdf_file(filename,'r')
        self.nSUB = ncIN.dimensions['nsub']
        self.nCOAST = ncIN.dimensions['ncoast']
        self.SUBLIST = ncIN.sublist[:-1].split(",")
        self.COASTLIST=ncIN.coastlist[:-1].split(",")
        ncIN.close()

    def read_validation_file(self,filename):
        ncIN = NC.netcdf_file(filename,'r')
        NUMBER = ncIN.variables['number'].data.copy()
        BIAS   = ncIN.variables['BGC_CLASS4_CHL_BIAS_SURF_BASIN'].data.copy()
        RMSE   = ncIN.variables['BGC_CLASS4_CHL_RMS_SURF_BASIN'].data.copy()
        MODEL  = ncIN.variables['MODEL_MEAN'].data.copy()
        SAT    = ncIN.variables['SAT___MEAN'].data.copy()
        BIASLOG= ncIN.variables['BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG'].data.copy()
        RMSELOG= ncIN.variables['BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG'].data.copy()     
        ncIN.close()
        return NUMBER,BIAS,RMSE, MODEL, SAT, BIASLOG, RMSELOG
    
    def readfiles(self):
        
        BIAS    =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        RMSE    =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        NUMBER  =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        MODEL   =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        SAT     =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        BIASLOG =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        RMSELOG =np.zeros((self.nFrames,self.nSUB,self.nCOAST))
        
        for iFrame, filename in enumerate(self.filelist):
            time = datetime.datetime.strptime(filename[-11:],'%Y%m%d.nc')
            self.timelist.append(time)
            number, bias, rmse, model, sat, logbias, logrmse= self.read_validation_file(filename)
            NUMBER[iFrame,:,:] = number
            BIAS[iFrame,:,:]   = bias
            RMSE[iFrame,:,:]   = rmse
            MODEL[iFrame,:,:]  = model
            SAT[iFrame,:,:]    = sat
            BIASLOG[iFrame,:,:] = logbias
            RMSELOG[iFrame,:,:] = logrmse
        self.bias = BIAS
        self.number=NUMBER
        self.rmse = RMSE
        self.model = MODEL
        self.sat   = SAT
        self.biaslog = BIASLOG
        self.rmselog = RMSELOG

    def append_dir(self,dirname):
        dirlist=glob.glob(dirname + "Validation_" + self.search_type + "*")
        if len(dirlist)> 0:
            self.append_file(dirlist[0])
        

    def append_file(self,filename):
        number, bias, rmse, model, sat, logbias, logrmse= self.read_validation_file(filename)
        time = datetime.datetime.strptime(filename[-11:],'%Y%m%d.nc')
        
        
        self.number = np.concatenate((self.number , number.reshape((1,self.nSUB, self.nCOAST))), axis=0)
        self.bias   = np.concatenate((self.bias   ,   bias.reshape((1,self.nSUB, self.nCOAST)) ), axis=0)
        self.rmse   = np.concatenate((self.rmse   ,   rmse.reshape((1,self.nSUB, self.nCOAST)) ), axis=0)
        self.model  = np.concatenate((self.model  ,  model.reshape((1,self.nSUB, self.nCOAST)) ), axis=0)
        self.sat    = np.concatenate((self.sat    ,    sat.reshape((1,self.nSUB, self.nCOAST)) ), axis=0)
        self.biaslog= np.concatenate((self.biaslog,logbias.reshape((1,self.nSUB, self.nCOAST)) ), axis=0)
        self.rmselog= np.concatenate((self.rmselog,logrmse.reshape((1,self.nSUB, self.nCOAST)) ), axis=0)
        self.timelist.append(time)
        self.nFrames = self.nFrames + 1
    
    
    def plotdata(self,VAR, sub,coast):
        '''
        VAR must be a 3D field of this class, such as bias, number, ...
        If the number of matchups is less than 50 - at 1/24 grid - the statistic won't be plotted
        '''

        isub = self.SUBLIST.index(sub)
        icoast= self.COASTLIST.index(coast)
        numb = self.number[:,isub, icoast]
        ii  = numb < 50
        y = VAR[:,isub,icoast]
        y[ii] = np.nan
        return self.timelist, y
