import netCDF4
import numpy as np
import datetime
import os
import matplotlib.pyplot as pl
from bitsea.commons.utils import addsep
from bitsea.instruments.instrument import Instrument, Profile
from scipy.optimize import curve_fit
from bitsea.instruments.var_conversions import FLOATVARS as conversion

mydtype= np.dtype([
          ('file_name','U200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','U17'),
          ('parameters','U200'),
          ('prof_flags' ,'U200')]
          )
GSS_DEFAULT_LOC ='/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C'
ONLINE_REPO = addsep(os.getenv("ONLINE_REPO",GSS_DEFAULT_LOC))
FloatIndexer=addsep(ONLINE_REPO) + "PPCON/Float_Index.txt"
INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)

nFiles=INDEX_FILE.size
for iFile in range(nFiles):
    INDEX_FILE['parameters'][iFile] = INDEX_FILE['parameters'][iFile].replace("NITRATE_PPCON","NITRATE").replace("CHLA_PPCON","CHLA").replace("BBP700_PPCON","BBP700")

class BioFloatProfile(Profile):
    def __init__(self, time, lon, lat, my_float, available_params , profile_flags ,mean=None):
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        if available_params[0]==" ":
           self.available_params = available_params[1:]
        else:
           self.available_params = available_params
        self.profile_flags= profile_flags
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

    def read(self,var,var_mod=None,sourcedata='best'):
        '''
        Reads profile data from file. Wrapper for BioFloat.read()

        Arguments:
        * var *  string
        * read_adjusted * IS NOT USED, but we leave it here because there is a lot of calls using it
                          Once all bit.sea code will use superfloat instead of lovbio_float, we'll remove it

        Returns 3 numpy arrays: Pres, Profile, Qc
        ProfileType codes= "best" take all profiles || "insitu" take insitu only || "ppcon" take ppcon only
        ProfileType default=> "best"
        '''

        return self._my_float.read(var,var_mod,sourcedata)


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

    def __init__(self,lon,lat,time,filename,available_params, profile_flags):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename
        if available_params[0]==" ":
            self.available_params = available_params[1:]
        else:
            self.available_params = available_params
        self.profile_flags    = profile_flags.replace(" ","")
        wmo, cycle = os.path.basename(filename).rsplit("_")
        self.wmo = wmo[2:]
        self.cycle = int(cycle[:3])

    def __eq__(self,other):
        if isinstance(other, BioFloat):
            if (self.filename  == other.filename):
                return (self.lon == other.lon ) & (self.lat == other.lat) & (self.time == other.time)
            else:
                return False
        else:
            return False



    def status_profile(self,var):
        """return code or extended nomenclature for insitu, ppcon, both
        """
        li = list(self.available_params.split(" "))
        lf = self.profile_flags
        iParam = li.index(var)
        flags_code = lf[iParam]
        return (flags_code)


    def has_insitu(self, var):
        if self.status_profile(var) in ['I','B']: return True
        return False
    def has_ppcon(self,  var):
        if self.status_profile(var) in['P','B']: return True
        return False
    def has_insitu_and_ppcon(self, var):
        if self.status_profile(var) == 'B': return True
        return False


    def read_raw(self , var , flag_code):
        '''
        Reads data from file
        Returns 3 numpy arrays: Pres, Profile, Qc
        '''
        ncIN = netCDF4.Dataset(self.filename,'r')
        if flag_code in ['I','B'] :
            Pres    = np.array(ncIN.variables['PRES_'+var])
            Profile = np.array(ncIN.variables[        var])
            Qc      = np.array(ncIN.variables[var    + "_QC"])
        else:
            Pres    = np.array(ncIN.variables['PRES_'+var+'_PPCON'])
            Profile = np.array(ncIN.variables[        var+'_PPCON'])
            Qc      = np.array(ncIN.variables[var+'_PPCON' + "_QC"])

        ncIN.close()
        return Pres, Profile, Qc


    def read(self,var,var_mod=None,sourcedata='best'):
        '''
        With sourcedata='best' gives priority to insitu

        '''

        assert sourcedata in ["best","insitu","ppcon"]

        if sourcedata=='best':
            flag_code = self.status_profile(var) # get flag insitu ppcon both none
            Pres, Profile, Qc = self.read_raw(var, flag_code)

        if sourcedata=="insitu":
            assert self.has_insitu(var)
            Pres, Profile, Qc  =  self.read_raw(var, 'I')
            if (Qc <0).any() : sys.exit("the sourcedata flag is uncorrect")

        if sourcedata=="ppcon":
            assert self.has_ppcon(var)
            Pres, Profile, Qc = self.read_raw(var, 'P')
            if (Qc >0).any() : sys.exit("the sourcedata flag is uncorrect")   

        if var_mod is None:
            return Pres, Profile, Qc

        #assert var_mod in ['P_c', 'POC']

        ii=(Pres >= 180) & (Pres <= 200)
        if (var_mod=='P_c'):
            bbp470 = Profile * ( 470.0/ 700)**(-0.78)# [m-1]
            Profile = 12128 * bbp470 + 0.59 # Conversion by Bellacicco 201?
            if ii.sum() > 0 :
                shift=Profile[ii].mean()
                print( var_mod  +" : adding a shift of " + str(shift))
                Profile = Profile - shift
                ii=Profile<=0
                Profile[ii] = 0.0

        if (var_mod == "POC"):
            POC = Profile *  52779.37 - 3.57 # Bellacicco 2019
            if ii.sum() > 0 :
                shift=POC[ii].mean()
                print( "POC: adding a shift of " + str(shift))
                Profile = POC - shift
                ii=Profile<=0
                Profile[ii] = 0.0

     
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

        ncIN=netCDF4.Dataset(self.filename,'r')
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
            timestr          = INDEX_FILE['time'][iFile]
            lon              = INDEX_FILE['lon' ][iFile]
            lat              = INDEX_FILE['lat' ][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile]
            available_params = INDEX_FILE['parameters'][iFile]
            profile_flags    = INDEX_FILE['prof_flags'][iFile]
            float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
            if filename.endswith(thefilename):
                return BioFloat(lon,lat,float_time,filename,available_params, profile_flags)
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
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        profile_flags    = INDEX_FILE['prof_flags'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
        filename = ONLINE_REPO + 'PPCON/' + filename

        if var is None :
            VarCondition = True
        else:
            VarCondition = var in available_params

        if VarCondition:
            if T.contains(float_time) and region.is_inside(lon, lat):
                thefloat = BioFloat(lon,lat,float_time,filename,available_params, profile_flags)
                selected.append(BioFloatProfile(float_time,lon,lat, thefloat,available_params,profile_flags,None))

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
    filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C/PPCON/5907088/SR5907088_016.nc"
    F=BioFloat.from_file(filename)
    print("profile_flags=",F.profile_flags)
    print("status_profile=", F.status_profile('NITRATE'))
    print("has_insitu", F.has_insitu('NITRATE'))
    F.read('DOXY')
    F.read('DOXY',sourcedata='insitu')
    Pres_P,Value_P, _ =F.read('CHLA',sourcedata='ppcon')
    Pres_I,Value_I, Qc=F.read('CHLA',sourcedata='insitu')
    Pres,Value, Qc=F.read('CHLA',sourcedata='best')
    assert (Value == Value_I).all()

    
    from bitsea.basins.region import Rectangle
    from bitsea.commons.time_interval import TimeInterval

    var = 'NITRATE'
    TI = TimeInterval('20240201','20240301','%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILE_LIST=FloatSelector(var, TI, R)
    for p in PROFILE_LIST:
        print(p._my_float.has_insitu('NITRATE'))
        Pres, Values, Qc = p.read('NITRATE')

    import sys
    sys.exit()
    
    Pres,V, Qc = F.read(var)
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
