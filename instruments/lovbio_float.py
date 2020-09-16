import scipy.io.netcdf as NC
import numpy as np
import datetime
import os
import matplotlib.pyplot as pl
import seawater as sw
from commons.utils import addsep
from var_conversions import LOVFLOATVARS as conversion

from instrument import Instrument, Profile
from mhelpers.pgmean import PLGaussianMean
meanObj = PLGaussianMean(5,1.0)

mydtype= np.dtype([
          ('file_name','S200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','S17'),
          ('parameters','S200')] )
GSS_DEFAULT_LOC = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/"
ONLINE_REPO = addsep(os.getenv("ONLINE_REPO",GSS_DEFAULT_LOC))
FloatIndexer=addsep(ONLINE_REPO) + "FLOAT_LOVBIO/Float_Index.txt"
is_default_V4C= ONLINE_REPO == GSS_DEFAULT_LOC

INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)


class BioFloatProfile(Profile):
    def __init__(self, time, lon, lat, my_float, available_params,mean=None):
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        self.available_params = available_params
        self.mean = mean
        self.has_adjusted = True

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

    def ID(self):
        return  self.name() + "_" + self.time.strftime("%Y%m%d_") + str(self.lon) + "_"+ str(self.lat)
    def reference_var(self,var):
        '''
        Returns the reference varname, for a given profile object and
        a ogstm model varname
        For BioFloats p.reference_var('O2o') returns 'DOXY'
        '''
        return conversion[var]

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
        '''
        If the variable is a char array, with fillvalue=''
        this function results in an array of '0'
        '''
        varObj = ncObj.variables[var]
        try:
            fillvalue = varObj._FillValue
        except:
            try:
                fillvalue = varObj.missing_value
            except:
                fillvalue = 99999
        M = varObj.data.copy()
        if M.dtype == np.dtype('S1'):
            M[M==fillvalue] = '0'
        else:
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
        PRES    = self.__merge_var_with_adjusted(ncIN, 'PRES')

        if ncIN.variables.has_key(var + "_ADJUSTED"):
            M_ADJ   = self.__fillnan(ncIN, var + "_ADJUSTED")
            QC      = self.__fillnan(ncIN, var +"_ADJUSTED_QC")
        else:
            M_ADJ   = np.zeros_like(PRES)*np.nan
            QC      = np.zeros_like(PRES)*np.nan
        M       = self.__fillnan(ncIN, var )
        #M      = self.__merge_var_with_adjusted(ncIN, var) # to have not only adjusted

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



        # Elimination of negative pressures or nans
        nanPres = np.isnan(rawPres)
        rawPres[nanPres] = 1 # just for not complaining
        badPres    = (rawPres<=0) | nanPres
        #goodQc     = (rawQc == '2' ) | (rawQc == '1' )
        nanProfile = (np.isnan(rawProfile))
        rawProfile [nanProfile] = 1  # just for not complaining
        badProfile = (rawProfile>1e+20) | nanProfile 
        bad = badPres | badProfile #| (~goodQc )


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
        if var=='CHLA': prof = prof*0.5
        if (var=='SR_NO3') & read_adjusted :
            #prof = prof * 0.98 + 0.6
            # New adjustement following Mignot et al. (2019)
            prof = prof * 1.04 + 0.46 
            ii = (prof < 0) & (pres < 50)
            prof[ii] = 0.05
            ii = prof > 0
            pres = pres[ii]
            prof = prof[ii]
            qc   =   qc[ii]
        if (var=='DOXY'):
            #prof = self.convert_oxygen(pres, prof)
            ii = (prof > 140) & (prof<280)
            pres = pres[ii]
            prof = prof[ii]
            qc   =   qc[ii]
            

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

    def convert_oxygen(self,doxypres,doxyprofile):
        '''
        from micromol/Kg to  mmol/m3
        '''
        if doxypres.size == 0: return doxyprofile
        Pres, temp, Qc = self.read_raw("TEMP",False)
        Pres, sali, Qc = self.read_raw("PSAL",False)
        density = sw.dens(sali,temp,Pres)
        density_on_zdoxy = np.interp(doxypres,Pres,density)
        return doxyprofile * 1000./density_on_zdoxy
        
        
        

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
        nFiles=INDEX_FILE.size
        for iFile in range(nFiles):
            timestr          = INDEX_FILE['time'][iFile]
            lon              = INDEX_FILE['lon' ][iFile]
            lat              = INDEX_FILE['lat' ][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile]
            available_params = INDEX_FILE['parameters'][iFile]
            float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
            if ONLINE_REPO + "FLOAT_LOVBIO/" + thefilename == filename :
                return BioFloat(lon,lat,float_time,filename,available_params)
        return None

def profile_gen(lon,lat,float_time,filename,available_params):
    if not is_default_V4C:
        filename = ONLINE_REPO + "FLOAT_LOVBIO/" + filename
    thefloat = BioFloat(lon,lat,float_time,filename,available_params)
    return BioFloatProfile(float_time,lon,lat, thefloat,available_params,meanObj)
def FloatSelector(var, T, region):
    '''
    Arguments:
       * var *    is a string indicating variable, 
                  if var is None, no selection is done about variable
       * T   *    is a TimeInterval instance or a timerequestors.Clim_season instance
                  or whatever object having a contains() method working as the TimeInteval one does.
       * region * is an instance of Region or its derived (Polygon, Basin, ...)
       
    Returns:
        a list of BioFloatProfile objects.
    Caveats:
       In order to work on dataset different from the cineca DRES archive
       /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/
       remember to define the environment variable ONLINE_REPO
       export ONLINE_REPO=/some/path/with/ COPERNICUS/  FLOAT_BIO/  FLOAT_LOVBIO/  SAT/
    '''

    nFiles=INDEX_FILE.size
    selected = []
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
        if not is_default_V4C :
            filename = ONLINE_REPO + "FLOAT_LOVBIO/" + filename

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
    wmo_set=set()
    for p in Profilelist:
        wmo_set.add(p._my_float.wmo)
    return list(wmo_set)

def filter_by_wmo(Profilelist,wmo):
    '''

    Subsetter, filtering by wmo

     Arguments:
      * Profilelist * list of Profile objects
      * wmo         * string

      Returns:
        a list of Profile Objects
    '''

    return [p for p in Profilelist if p._my_float.wmo == wmo]

def remove_bad_sensors(Profilelist,var):
    '''

    Subsetter, filtering out bad sensors for that var

     Arguments:
      * Profilelist * list of Profile objects
      * var         * string

      Returns:
        a list of Profile Objects
    '''
 
    OUT_N3n = ["6903197","6901767","6901773","6901771"]
    OUT_O2o = ["6901766","6901510"]

    if ( var == 'SR_NO3' ):
        return [p for p in Profilelist if p.name() not in OUT_N3n]

    if ( var == 'DOXY' ):
        return [p for p in Profilelist if p.name() not in OUT_O2o]

    return Profilelist

def from_coriolis_profile(coriolis_profile, verbose=True):
    '''
    Arguments:
    * lov_profile * a lov profile objects
    * verbose     * logical, used to print lov files that don't have a corresponding in coriolis

    Returns:
    *  p * a LOV BioFloatProfile object
    '''
    wmo = coriolis_profile._my_float.wmo
    INDEXES=[]
    for iFile, filename in enumerate(INDEX_FILE['file_name']):
        if filename.startswith(wmo):
            INDEXES.append(iFile)
    A = INDEX_FILE[INDEXES]
    nFiles = len(A)
    DELTA_TIMES = np.zeros((nFiles,), np.float32)
    for k in range(nFiles):
        float_time =datetime.datetime.strptime(A['time'][k],'%Y%m%d-%H:%M:%S')
        deltat = coriolis_profile.time - float_time
        DELTA_TIMES[k] = deltat.total_seconds()
    min_DeltaT = np.abs(DELTA_TIMES).min()
    if min_DeltaT > 3600*3 :
        if verbose: print "no LOV file corresponding to "  + coriolis_profile._my_float.filename
        return None
    F = (A['lon'] - coriolis_profile.lon)**2 + (A['lat'] - coriolis_profile.lat)**2 +  DELTA_TIMES**2
    iFile = F.argmin()
    timestr          = A['time'][iFile]
    lon              = A['lon' ][iFile]
    lat              = A['lat' ][iFile]
    filename         = A['file_name'][iFile]
    available_params = A['parameters'][iFile]
    float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
    return profile_gen(lon, lat, float_time, filename, available_params)



if __name__ == '__main__':
    from basins.region import Rectangle
    from commons.time_interval import TimeInterval

    var = 'DOXY'
#    var = 'SR_NO3'
    TI = TimeInterval('20140120','20170130','%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILE_LIST=FloatSelector(var, TI, R)
    
    print len(PROFILE_LIST)
    print len(remove_bad_sensors(PROFILE_LIST,var))



    for ip, p in enumerate(PROFILE_LIST):
        continue
        F = p._my_float
        Pres,V, Qc = F.read(var, read_adjusted=False)
        ii =~np.isnan(V)
        V = V[ii]
        if len(V)> 0:
            print V.max(), V.min()
            if V.min() < -976:
                break


    for p in PROFILE_LIST[:1]:
        PN,N, Qc = p.read(var,read_adjusted=True)
        TheFloat = p._my_float
        PN,N,Qc = TheFloat.read(var,    True)
        PS,S,Qc = TheFloat.read('PSAL', True)
        PT,T,Qc = TheFloat.read('TEMP', True)

    wmo_list= get_wmo_list(PROFILE_LIST)
    for wmo in wmo_list:
        sublist = filter_by_wmo(PROFILE_LIST, wmo)
