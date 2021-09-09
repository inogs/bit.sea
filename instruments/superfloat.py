import scipy.io.netcdf as NC
import numpy as np
import datetime
import os
import matplotlib.pyplot as pl
from commons.utils import addsep
from instruments.instrument import Instrument, Profile
from scipy.optimize import curve_fit
from instruments.var_conversions import FLOATVARS as conversion

mydtype= np.dtype([
          ('file_name','S200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','S17'),
          ('parameters','S200')] )
GSS_DEFAULT_LOC = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V7C/"
ONLINE_REPO = addsep(os.getenv("ONLINE_REPO",GSS_DEFAULT_LOC))
FloatIndexer=addsep(ONLINE_REPO) + "SUPERFLOAT/Float_Index.txt"
INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)

class BioFloatProfile(Profile):
    def __init__(self, time, lon, lat, my_float, available_params,mean=None):
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        self.available_params = available_params
        self.mean = mean
        self.has_adjusted=False

    def __eq__(self, other):
        if isinstance(other, BioFloatProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                return self._my_float == other._my_float
            else:
                return False
        else:
            return False

    def read(self,var,var_mod=None):
        '''
        Reads profile data from file. Wrapper for BioFloat.read()

        Arguments:
        * var *  string
        * read_adjusted * IS NOT USED, but we leave it here because there is a lot of calls using it
                          Once all bit.sea code will use superfloat instead of lovbio_float, we'll remove it

        Returns 3 numpy arrays: Pres, Profile, Qc '''

        return self._my_float.read(var,var_mod)


    def read_fitted(self,var, func):
        '''
        This is used to calculate the radiometric profile values at the surface.
        Use an exponential fit and extrapolate values to the surface.
        Returns 3 numpy arrays: Pres, Profile, Qc with values extrapolated to 0 m
        '''

        return self._my_float.read_fitted(var, func, mean=self.mean)


    def name(self):
        '''returns a string, the wmo of the associated BioFloat.
        '''
        return self._my_float.wmo

    def ID(self):
        return  self.name() + "_" + self.time.strftime("%Y%m%d-%H:%M:%S_") + str(self.lon) + "_"+ str(self.lat)
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
        #istart=filename.index("/",filename.index('SUPERFLOAT/'))
        #iend  =filename.index("/",istart+1)
        wmo, cycle = os.path.basename(filename).rsplit("_")
        self.wmo = wmo[2:]
        self.cycle = int(cycle[:3])
        #self.wmo = filename[istart+1:iend]
        #cycle = os.path.splitext(os.path.basename(filename))[0].rsplit("_")[1]
        #self.cycle = int(cycle)


    def __eq__(self,other):
        if isinstance(other, BioFloat):
            if (self.filename  == other.filename):
                return (self.lon == other.lon ) & (self.lat == other.lat) & (self.time == other.time)
            else:
                return False
        else:
            return False



    def read_raw(self,var):
        '''
        Reads data from file
        Returns 3 numpy arrays: Pres, Profile, Qc
        '''
        ncIN=NC.netcdf_file(self.filename,'r')
        Pres    = ncIN.variables['PRES_'+var].data.copy()
        Profile = ncIN.variables[        var].data.copy()
        Qc      = ncIN.variables[var + "_QC"].data.copy()
        ncIN.close()
        return Pres, Profile, Qc


    def read(self,var,var_mod=None):

        Pres, Profile, Qc = self.read_raw(var)

        if var_mod is None:
#           print "var_mod is " + np.str(var_mod) 
           return Pres, Profile, Qc

        ii=(Pres >= 400) & (Pres <= 500) 
        if (var_mod=='P_c'):
            bbp470 = Profile * ( 470.0/ 700)**0.78# [m-1]
            Profile = 12128 * bbp470 + 0.59 # Conversion by Bellacicco 201?
            shift=Profile[ii].mean()
            Profile = Profile - shift
            ii=Profile<=0
            Profile[ii] = 0.0

#        if (var_mod == "POC"): 

     
        return Pres, Profile, Qc


    def read_fitted(self, var, func, mean=None):
        '''
        This is used to calculate the radiometric profile values at the surface.
        Use an exponential fit and extrapolate values to the surface.
        Returns the irradiance values at the surface.

        Reads profile data from file, applies a rarefaction and optionally a filter to the data

        Returns the measurement's surface value

        '''

        raw_pres, raw_prof   = self.read_raw(var)
        pres,index           = np.unique(raw_pres,return_index=True)
        prof                 = raw_prof[index]

        if pres[0] < 1.5 and pres[3] < 10.: 
        
            popt,_   = curve_fit(func, pres, prof)
            pres_new = np.insert(pres, 0, 0.)  if not pres[0] == 0. else pres # Add z=0 m
            prof_new = func(pres_new, *popt)

        else:
            pres_new = pres
            prof_new = prof

        qc       = np.ones_like(pres_new)*2

        #import matplotlib.pyplot as plt
        #plt.plot(prof, -pres, 'm', label='Profile')
        #plt.plot(prof_new, -pres_new, 'b', label='Curve fit')
        #plt.legend(loc='best')
        #plt.show()
        
        if pres_new.size ==0:
            return pres_new, prof_new, qc

        if mean == None:
            if BioFloat.default_mean != None:
                return pres_new, BioFloat.default_mean.compute(prof_new, pres_new), qc
            else:
                return pres_new, prof_new, qc
        else:
            return pres_new, mean.compute(prof_new, pres_new), mean.compute(qc,pres_new)


    def origin(self,var):
        '''
        Arguments:
        * var * string
        Returns:
        a tuple of two strings:
           * origin      * string, "lov" or "Coriolis"
           * file_origin * string , path of the source file
        '''
        class Info():
            def __init__(self):
                self.file_orig   = "None"
                self.correction = "None"
                self.shift      = "None"
                self.status_var = 'n'
            def __repr__(self):
                return "status_var = %s" %(self.status_var) 

        info = Info()

        ncIN=NC.netcdf_file(self.filename,'r')
        info.status_var=ncIN.variables[var].status_var
        info.file_orig=ncIN.file_origin
        
        if var=='NITRATE':
            info.shift=ncIN.variables[var].bottomshift
            info.correction=ncIN.variables[var].correction
        ncIN.close()
        
        return info


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
            timestr          = INDEX_FILE['time'][iFile].decode()
            lon              = INDEX_FILE['lon' ][iFile]
            lat              = INDEX_FILE['lat' ][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile].decode()
            available_params = INDEX_FILE['parameters'][iFile].decode()
            float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
            if filename.endswith(thefilename):
                return BioFloat(lon,lat,float_time,filename,available_params)
        return None


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
       export ONLINE_REPO=/some/path/with/ COPERNICUS/  FLOAT_BIO/  FLOAT_LOVBIO/  SAT/ SUPERFLOAT/
    '''

    nFiles=INDEX_FILE.size
    selected = []
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile].decode()
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile].decode()
        available_params = INDEX_FILE['parameters'][iFile].decode()
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
        filename = ONLINE_REPO + "SUPERFLOAT/" + filename

        if var is None :
            VarCondition = True
        else:
            VarCondition = var in available_params

        if VarCondition:
            if T.contains(float_time) and region.is_inside(lon, lat):
                thefloat = BioFloat(lon,lat,float_time,filename,available_params)
                selected.append(BioFloatProfile(float_time,lon,lat, thefloat,available_params,None))

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

    Subsetter, filtering out bad sensors for that var.
    At the moment that method does not do anything, because bad sensors have been already removed
    by superfloat generators.

    In case of new expert evaluation, if we decide to remove another bad sensor,
    we can add here the part of code to have quickly the desired result.
    Then, the procedure to avoid bad sensors in superfloat dataset is
    - adding bad float in superfloat_${var}.py remove_bad_sensor()
    - remove files from SUPERFLOAT DATASET
    - launch dump_index.py without Float_Indexer.txt input


     Arguments:
      * Profilelist * list of Profile objects
      * var         * string

      Returns:
        a list of Profile Objects
    '''

    return Profilelist


if __name__ == '__main__':
    from basins.region import Rectangle
    from commons.time_interval import TimeInterval

    var = 'TEMP'
    TI = TimeInterval('20120101','20170130','%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILE_LIST=FloatSelector(var, TI, R)
    filename="/gpfs/scratch/userexternal/gbolzon0/V7C/DEBUG_SUPERFLOAT/ONLINE/SUPERFLOAT/6902969/SR6902969_001.nc"
    F=BioFloat.from_file(filename)
    import sys
    sys.exit()

#    print len(PROFILE_LIST)

    for p in PROFILE_LIST[100:200]:
        Pres,V, Qc = p.read(var)
#        if Pres.min()>0:
#            print Pres.min()

    wmo_list= get_wmo_list(PROFILE_LIST)
    for wmo in wmo_list:
        sublist = filter_by_wmo(PROFILE_LIST, wmo)
