import scipy.io.netcdf as NC
import numpy as np
import datetime
import os
import pylab as pl

from instrument import Instrument, Profile
from mhelpers.pgmean import PLGaussianMean
meanObj = PLGaussianMean(5,1.0)



class BioFloatProfile(Profile):
    def __init__(self, time, lon, lat, my_float, available_params,mean=None):
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        self.available_params = available_params
        self.mean = mean

    def __eq__(self, other):
        if isinstance(other, BioFloatProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                return self._my_float == other._my_float
            else:
                return False
        else:
            return False

    def read(self,var,read_adjusted):
        '''
        Reads profile data from file. Wrapper for BioFloat.read()

        Takes var as string
              read_adjusted as logical
        Returns 3 numpy arrays: Pres, Profile, Qc '''

        return self._my_float.read(var, mean=self.mean,read_adjusted=read_adjusted)

    def name(self):
        '''returns a string, the wmo of the associated BioFloat.
        '''
        return self._my_float.wmo



class BioFloat(Instrument):

    default_mean = None

    def __init__(self,lon,lat,time,filename,available_params):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename
        self.available_params = available_params
        istart=filename.index("/",filename.index('FLOAT_LOVBIO'))
        iend  =filename.index("/",istart+1)
        self.wmo = filename[istart+1:iend]
        cycle = os.path.basename(filename).rsplit("_")[2]
        self.cycle = int(cycle)

    def __eq__(self,other):
        if isinstance(other, BioFloat):
            if (self.filename  == other.filename):
                return (self.lon == other.lon ) & (self.lat == other.lat) & (self.time == other.time)
            else:
                return False
        else:
            return False

    def __searchVariable_on_parameters(self, var):
        '''
        returns index of profile on which variable has to be read,
        -1 if fails (PARAMETERS variable is blank)

        '''

        ncIN = NC.netcdf_file(self.filename,'r')
        N_PROF    = ncIN.dimensions['N_PROF']
        N_PARAM   = ncIN.dimensions['N_PARAM']
        PARAMETER = ncIN.variables['PARAMETER'].data.copy()
        ncIN.close()

        if PARAMETER.tostring().rstrip() == '':
            return -1
        else:
            for iprof in range(N_PROF):
                for iparam in range(N_PARAM):
                    s = PARAMETER[iprof,0,iparam,:].tostring().rstrip()
                    if s==var:
                        return iprof

    def __fillnan(self, ncObj,var):
        varObj = ncObj.variables[var]
        fillvalue = varObj._FillValue
        M = varObj.data.copy()
        M[M==fillvalue] = np.NaN;
        return M

    def __merge_profile_with_adjusted(self,profile, profile_adj):
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

    def __merge_var_with_adjusted(self, ncObj,var):
        N_PROF= ncObj.dimensions['N_PROF']
        M     = self.__fillnan(ncObj, var)
        M_ADJ = self.__fillnan(ncObj, var + "_ADJUSTED")
        M_RES = M
        for iprof in range(N_PROF):
            M_RES[iprof,:] = self.__merge_profile_with_adjusted(M[iprof,:], M_ADJ[iprof,:])

        return M_RES



    def read_very_raw(self,var):
        '''
        Reads data from file
        Returns 4 numpy arrays: Pres, Profile, Profile_adjusted, Qc
        '''
        
        iProf = 0#self.__searchVariable_on_parameters(var); #print iProf
        ncIN=NC.netcdf_file(self.filename,'r')

        if iProf== -1 :
            M = self.__merge_var_with_adjusted(ncIN, var)
            N_PROF= ncIN.dimensions['N_PROF']
            N_MEAS = np.zeros((N_PROF),np.int32)
            for iprof in range(N_PROF):
                N_MEAS[iprof] = (~np.isnan(M[iprof,:])).sum()
            #per il momento cerco il massimo
            iProf = N_MEAS.argmax()
            print "iprof new value", iProf


        M_ADJ   = self.__fillnan(ncIN, var + "_ADJUSTED")
        M       = self.__fillnan(ncIN, var )
        #M      = self.__merge_var_with_adjusted(ncIN, var) # to have not only adjusted
        PRES    = self.__merge_var_with_adjusted(ncIN, 'PRES')
        QC      = self.__fillnan(ncIN, var +"_ADJUSTED_QC")
        Profile     =     M[iProf,:]
        Profile_adj = M_ADJ[iProf,:]
        Pres        =  PRES[iProf,:]
        Qc          =    QC[iProf,:]
        ncIN.close()
        
        return Pres, Profile,Profile_adj, Qc

    def read_raw(self,var,read_adjusted=True):
        '''
        Reads a clean profile data
         - without nans
         - without repetition of the same pressure
        Inputs:
          var (string)
          read_adjusted (logical)
        Returns 3 numpy arrays: Pres, Profile, Qc
        '''
        rawPres, rawProfile, rawProfile_adj, rawQc = self.read_very_raw(var)

        
        if read_adjusted:
            rawProfile = rawProfile_adj
        else:
            for i in range(len(rawQc)): rawQc[i]='2' # to force goodQc = True


        # Elimination of negative pressures or nans
        nanPres = np.isnan(rawPres)
        rawPres[nanPres] = 1 # just for not complaining
        badPres    = (rawPres<=0) | nanPres
        goodQc     = (rawQc == '2' ) | (rawQc == '1' )
        badProfile = np.isnan(rawProfile)
        bad = badPres | badProfile | (~goodQc )


        Pres    =    rawPres[~bad]
        Profile = rawProfile[~bad]
        Qc      =     (rawQc[~bad]).astype(np.int)

        uniquePres,index=np.unique(Pres,return_index=True)
        uniqueProfile =  Profile[index]
        uniqueQc      =       Qc[index]

        return uniquePres, uniqueProfile, uniqueQc

    def rarefy(self,Pres,minimumdelta):
        '''
        Used to reduce the amount of data 
        ( e.g the profiles having data every 0.2m ) to match with model profiles.
        minimumdelta is expressed in meters.

        Returns a numpy array of indexes that ensure that
        Pres[indexes] has steps greater then minimumdelta

        '''
        indexes=[0]
        lastind=0
        for i,p in enumerate(Pres):
            if p >= Pres[lastind]+minimumdelta:
                indexes.append(i)
                lastind = i

        return np.array(indexes)

    def read(self, var, mean=None, read_adjusted=True):
        '''

        Reads profile data from file, applies a rarefaction and optionally a filter to the data

        Takes var as string
              read_adjusted as logical
        Returns 3 numpy arrays: Pres, Profile, Qc
        '''
        pres, prof, qc = self.read_raw(var,read_adjusted)
        if pres.size ==0:
            return pres, prof, qc

        ii = self.rarefy(pres, 2.0)
        pres = pres[ii]
        prof = prof[ii]
        qc   =   qc[ii]

        if mean == None:
            if BioFloat.default_mean != None:
                return pres, BioFloat.default_mean.compute(prof, pres), qc
            else:
                return pres, prof, qc
        else:
            return pres, mean.compute(prof, pres), mean.compute(qc,pres)

    def basicplot(self,Pres,profile):
        pl.figure()
        pl.plot(profile,Pres)
        pl.gca().invert_yaxis()
        pl.show(block=False)

    def plot(self,Pres,profile,fig=None,ax=None,**kwargs):
        '''
    Args:
        - * Pres    * : array of pressure
        - * profile * : array of profile values
        - * fig     * : a reference to a Figure object, if None a new Figure will be created.
        - * ax      * : a reference to an Axes object, if None a new axis will be created.

    Returns :
        A figure and an axis object

    Examples:
        fig, ax = f.plot(pres,profile)
        fig, ax = f.plot(pres,profile,fig,ax)
        fig, ax = f.plot(pres,profile,fig,ax, linestyle = 'dashed', linewidth = 2, color='green')

        '''
        if (fig is None) or (ax is None):
            fig , ax = pl.subplots()
        ax.plot(profile,Pres, **kwargs)
        if not ax.yaxis_inverted(): ax.invert_yaxis()
        return fig,ax

    def profiles(self, var, mean=None):
        return [BioFloatProfile(var, self.time, self.lon, self.lat, self, mean)]

    @staticmethod
    def from_file(filename):
        '''
        Returns the single Bio_Float instance corresponding to filename
        '''
        mydtype= np.dtype([
                  ('file_name','S200'),
                  ('lat',np.float32),
                  ('lon',np.float32),
                  ('time','S17'),
                  ('parameters','S200')] )

        FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_LOVBIO/Float_Index.txt"
        INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
        nFiles=INDEX_FILE.size
        for iFile in range(nFiles):
            timestr          = INDEX_FILE['time'][iFile]
            lon              = INDEX_FILE['lon' ][iFile]
            lat              = INDEX_FILE['lat' ][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile]
            available_params = INDEX_FILE['parameters'][iFile]
            float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
            if filename == thefilename :
                return BioFloat(lon,lat,float_time,filename,available_params)
        return None


def FloatSelector(var, T, region):
    '''
    Arguments:
       * var *    is a string indicating variable, 
                  if var is None, no selection is done about variable
       * T   *    is a TimeInterval instance
       * region * is an instance of Region or its derived (Polygon, Basin, ...)
       
    Returns:
        a list of BioFloatProfile objects.
    '''

    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )

    FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_LOVBIO/Float_Index.txt"

    INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
    nFiles=INDEX_FILE.size
    selected = []
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
            if T.contains(float_time) and region.is_inside(lon, lat):
                thefloat = BioFloat(lon,lat,float_time,filename,available_params)
                selected.append(BioFloatProfile(float_time,lon,lat, thefloat,available_params,meanObj))

    return selected

def get_wmo_list(Profilelist):
    '''
     Argument:
      * Profilelist * list of Profile objects

      Returns:
         a list of wmo strings
    '''
    raise NotImplementedError

def filter_by_wmo(Profilelist,wmo):
    '''

    Subsetter, filtering by wmo

     Arguments:
      * Profilelist * list of Profile objects
      * wmo         * string

      Returns:
        a list of Profile Objects
    '''
    raise NotImplementedError



if __name__ == '__main__':
    from basins.region import Rectangle
    from commons.time_interval import TimeInterval

    var = 'NITRATE'
    TI = TimeInterval('20150520','20150830','%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILE_LIST=FloatSelector(var, TI, R)

    for p in PROFILE_LIST[:1]:
        PN,N, Qc = p.read(var,read_adjusted=True)
        TheFloat = p._my_float
        PN,N,Qc = TheFloat.read(var,    True)
        PS,S,Qc = TheFloat.read('PSAL', True)
        PT,T,Qc = TheFloat.read('TEMP', True)

    wmo_list= get_wmo_list(PROFILE_LIST)
    for wmo in wmo_list:
        sublist = filter_by_wmo(PROFILE_LIST, wmo)
