import scipy.io.netcdf as NC
import numpy as np
import datetime
import os
import matplotlib.pyplot as pl
from bitsea.commons.utils import addsep
from var_conversions import FLOAT_OPT_VARS as conversion

from instrument import Instrument, Profile
from bitsea.mhelpers.pgmean import PLGaussianMean
meanObj = PLGaussianMean(5,1.0)



class BioFloatProfile(Profile):
    def __init__(self, time, lon, lat, my_float, available_params,mean=None):
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        self.available_params = available_params
        self.mean = mean
        self.has_adjusted = False

    def __eq__(self, other):
        if isinstance(other, BioFloatProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                return self._my_float == other._my_float
            else:
                return False
        else:
            return False

    def read(self,var):
        '''
        Reads profile data from file. Wrapper for BioFloat.read()

        Takes var as string
              read_adjusted as logical
        Returns 3 numpy arrays: Pres, Profile, Qc '''

        return self._my_float.read(var, mean=self.mean)

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
        istart=filename.index("/",filename.index('Float_OPT'))
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



    def read_raw(self,var):
        '''
        Reads data from file
        Returns 2 numpy arrays: Pres, Profile
        '''
        ncIN=NC.netcdf_file(self.filename,'r')
        Pres    = ncIN.variables['PRES_'+var].data.copy()
        Profile = ncIN.variables[        var].data.copy()
        ncIN.close()

        return Pres, Profile

    def read(self, var, mean=None):
        '''

        Reads profile data from file, applies a rarefaction and optionally a filter to the data

        Takes var as string
              read_adjusted as logical
        Returns 3 numpy arrays: Pres, Profile, Qc
        '''


        raw_pres, raw_prof   = self.read_raw(var)
        pres,index=np.unique(raw_pres,return_index=True)
        prof =  raw_prof[index]
        qc = np.ones_like(pres)*2
        #if var=='CHLA': prof = prof*0.5


        if pres.size ==0:
            return pres, prof, qc

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

        FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT/Float_Index.0.txt"
        INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
        nFiles=INDEX_FILE.size
        for iFile in range(nFiles):
            timestr          = INDEX_FILE['time'][iFile]
            lon              = INDEX_FILE['lon' ][iFile]
            lat              = INDEX_FILE['lat' ][iFile]
            thefilename      = INDEX_FILE['file_name'][iFile]
            available_params = INDEX_FILE['parameters'][iFile]
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
       export ONLINE_REPO=/some/path/with/ COPERNICUS/  FLOAT_BIO/  FLOAT_LOVBIO/  SAT/
    '''

    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )
    GSS_DEFAULT_LOC = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/"
    ONLINE_REPO = addsep(os.getenv("ONLINE_REPO",GSS_DEFAULT_LOC))
    FloatIndexer=addsep(ONLINE_REPO) + "Float_OPT/Float_Index.0.txt"


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
        filename = ONLINE_REPO + "Float_OPT/" + filename

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


if __name__ == '__main__':
    from bitsea.basins.region import Rectangle
    from bitsea.commons.time_interval import TimeInterval

    var = 'TEMP'
    TI = TimeInterval('20120101','20170130','%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILE_LIST=FloatSelector(var, TI, R)
    filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT/6901483/output_lovbio039b_040_00.nc"
    F=BioFloat.from_file(filename)
    import sys
    sys.exit()

    print len(PROFILE_LIST)

    for p in PROFILE_LIST:
        Pres,V, Qc = p.read(var)
        if Pres.min()>0:
            print Pres.min()

    wmo_list= get_wmo_list(PROFILE_LIST)
    for wmo in wmo_list:
        sublist = filter_by_wmo(PROFILE_LIST, wmo)
