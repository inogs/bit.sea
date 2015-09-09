import numpy as np
import datetime
import scipy.io.netcdf as NC
import pylab as pl
import os


from region import Region, Rectangle

class Time_Interval():
    def __init__(self,starttime="19500101", endtime="21000101", dateformat='%Y%m%d'):
        self.starttime=datetime.datetime.strptime(starttime,dateformat)
        self.end__time=datetime.datetime.strptime(endtime  ,dateformat)


class Bio_Float():
    def __init__(self,lon,lat,time,filename,available_params):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename
        self.available_params = available_params
        wmo,cycle=os.path.basename(filename).rsplit("_")
        self.wmo = wmo[2:]
        self.cycle = int(cycle[:3])
        
        self.T = 0
        self.S = 0

    def __searchVariable_on_parameters__(self,var):
        '''
        returns index of profile on which variable has to be read,
        -1 if fails (PARAMETERS variable is blank)

        '''

        ncIN = NC.netcdf_file(self.filename,'r')
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

    def __fillnan__(self, ncObj,var):
        varObj = ncObj.variables[var]
        fillvalue  =varObj._FillValue
        M     = varObj.data.copy()
        M[M==fillvalue] = np.NaN;
        return M

    def __merge_profile_with_adjusted__(self,profile, profile_adj):
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

    def __merge_var_with_adjusted__(self, ncObj,var):
        N_PROF= ncObj.dimensions['N_PROF']
        M     = self.__fillnan__(ncObj, var)
        M_ADJ = self.__fillnan__(ncObj, var + "_ADJUSTED")
        M_RES = M
        for iprof in range(N_PROF):
            M_RES[iprof,:] = self.__merge_profile_with_adjusted__(M[iprof,:], M_ADJ[iprof,:])

        return M_RES

    def read(self,var):
        '''
        Reads data from file
        '''
        iProf = self.__searchVariable_on_parameters__(var); #print iProf
        ncIN=NC.netcdf_file(self.filename,'r')

        if iProf== -1 :
            M = self.__merge_var_with_adjusted__(ncIN, var)
            N_PROF= ncIN.dimensions['N_PROF']
            N_MEAS = np.zeros((N_PROF),np.int32)
            for iprof in range(N_PROF):
                N_MEAS[iprof] = (~np.isnan(M[iprof,:])).sum()
            #per il momento cerco il massimo
            iProf = N_MEAS.argmax()
            print "iprof new value", iProf

        M = self.__merge_var_with_adjusted__(ncIN, var)
        PRES    = self.__merge_var_with_adjusted__(ncIN, 'PRES')
        Profile =  M[iProf,:]
        Pres    =  PRES[iProf,:]
        ncIN.close()

        # Elimination of negative pressures or nans
        badPres    = (Pres<=0) | (np.isnan(Pres))
        badProfile = np.isnan(Profile)
        bad = badPres | badProfile

        return Pres[~bad], Profile[~bad]
        #return Pres, Profile


    def plot(self,Pres,profile):
        pl.figure()
        pl.plot(profile,Pres)
        pl.gca().invert_yaxis()
        pl.show(block=False)

def FloatSelector(var, T, region):
    '''
    Arguments:
       var is a string indicating variable, 
          if var is None, no selection is done about variable
       T is as Time_Interval istance
       region is a region istance
    '''
    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )
    
    FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Index.txt"
    #FloatIndexer="MyFloat_Index.txt"
    INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
    nFiles=INDEX_FILE.size    
    SELECTED=[]
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
        if var is None :
            VarCondition = True
        else:
            VarCondition = var in available_params


        if VarCondition:
            if (float_time >= T.starttime) & (float_time <T.end__time):
                if region.is_inside(lon, lat):
                    SELECTED.append(Bio_Float(lon,lat,float_time,filename))
                      
    return SELECTED

if __name__ == '__main__':
    var = 'NITRATE'
    TI = Time_Interval('20150520','20150830','%Y%m%d')
    R = Rectangle(-6,36,30,46)
    
    FLOAT_LIST=FloatSelector(var, TI, R)

    for TheFloat in FLOAT_LIST[:1]:
        PN,N = TheFloat.read(var)
        PS,S = TheFloat.read('PSAL')
        PT,T = TheFloat.read('TEMP')

  
            
