import netCDF4
import numpy as np
import datetime
import os
from bitsea.commons.utils import addsep
from bitsea.instruments.var_conversions import FLOATVARS as conversion
import matplotlib
import matplotlib.pyplot as pl
from pathlib import Path
from bitsea.instruments.instrument import Instrument, Profile

CORIOLIS_DIR="CORIOLIS/"
mydtype= np.dtype([
          ('file_name','U200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','U17'),
          ('parameters','U200'),
          ('parameter_data_mode','U100')] )

GSS_DEFAULT_LOC = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C/"
ONLINE_REPO = addsep(os.getenv("ONLINE_REPO",GSS_DEFAULT_LOC))
FloatIndexer=addsep(ONLINE_REPO) + CORIOLIS_DIR + "Float_Index.txt"

INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)

class BioFloatProfile(Profile):
    def __init__(self, time, lon, lat, my_float, available_params):
        self.time = time
        self.lon = lon
        self.lat = lat
        self._my_float = my_float
        self.available_params = available_params
        self.has_adjusted = True

    def __eq__(self, other):
        if isinstance(other, BioFloatProfile):
            if (self.lon == other.lon) & (self.lat == other.lat) & (self.time == other.time):
                return self._my_float == other._my_float
            else:
                return False
        else:
            return False

    def read(self,var,read_adjusted=True):
        '''
        Reads profile data from file. Wrapper for BioFloat.read()

        Takes var as string
              read_adjusted as logical
        Returns 3 numpy arrays: Pres, Profile, Qc '''

        return self._my_float.read(var, read_adjusted=read_adjusted)

    def name(self):
        '''returns a string, the wmo of the associated BioFloat.
        '''
        return self._my_float.wmo

    def ID(self):
        return self.name() + "_" + self.time.strftime("%Y%m%d_") + str(self.lon) + "_"+ str(self.lat)
    def reference_var(self,var):
        '''
        Returns the reference varname, for a given profile object and
        a ogstm model varname
        For BioFloats p.reference_var('O2o') returns 'DOXY'
        '''
        return conversion[var]


class BioFloat(Instrument):

    def __init__(self,lon,lat,time,filename,available_params,parameter_data_mode):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename
        self.available_params = available_params
        wmo, cycle = os.path.basename(filename).rsplit("_")
        self.wmo = wmo[2:]
        self.cycle = int(cycle[:3])
        self.DATA_MODE = None
        self.PARAMETER_DATA_MODE = parameter_data_mode
        self.PARAMETER = None
        self.N_PARAM   = None

    def __eq__(self,other):
        if isinstance(other, BioFloat):
            if (self.filename  == other.filename):
                return (self.lon == other.lon ) & (self.lat == other.lat) & (self.time == other.time)
            else:
                return False
        else:
            return False
    def load_basics(self):

        ncIN = netCDF4.Dataset(self.filename,'r')
        self.N_PROF    = ncIN.dimensions['N_PROF'].size
        self.N_PARAM   = ncIN.dimensions['N_PARAM'].size
        self.PARAMETER = np.array(ncIN.variables['PARAMETER'])

        ncIN.close()
        _,_,nParams,_=self.PARAMETER.shape
        self.PARAMETERS=[]
        for iParam in range(nParams):
            Param_name=self.PARAMETER[0,0,iParam,:].tobytes().decode().replace(' ',"")
            self.PARAMETERS.append(Param_name)

    def _searchVariable_on_parameters(self, var):
        '''
        returns index of profile on which variable has to be read,
        -1 if fails (PARAMETERS variable is blank)

        '''
        if self.PARAMETER is None:
            self.load_basics()
        return self.PARAMETERS.index(var)


    def __fillnan(self, ncObj,var):
        varObj = ncObj.variables[var]
        fillvalue = varObj._FillValue
        M = np.array(varObj)
        if varObj.datatype=='S1':
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
        N_PROF= ncObj.dimensions['N_PROF'].size
        M     = self.__fillnan(ncObj, var)
        M_ADJ = self.__fillnan(ncObj, var + "_ADJUSTED")
        M_RES = M
        for iprof in range(N_PROF):
            M_RES[iprof,:] = self.__merge_profile_with_adjusted(M[iprof,:], M_ADJ[iprof,:])

        return M_RES

    def status_var(self,var):
        '''
        Argument: var, string
        Returns: string, 'R', 'A' or 'D' (real time or adjusted, or delay mode)
        '''
        param_list=self.available_params.rsplit(" ")
        iParam = param_list.index(var)
        return self.PARAMETER_DATA_MODE[iParam]
    def adjusted(self,var):
        if self.status_var(var) in ['A','D'] : return True
        return False
#         if self.PARAMETER is None:
#             self.load_basics()
#         iParam = self._searchVariable_on_parameters(var);
#         return self.PARAMETER_DATA_MODE[0,iParam]

    def read_very_raw(self,var):
        '''
        Reads data from file
        Returns 5 numpy arrays: Pres, Profile, Profile_adjusted, Qc, Qc_adjusted
        '''
        iProf = 0 #self._searchVariable_on_parameters(var)
        ncIN=netCDF4.Dataset(self.filename,'r')
        PRES    = self.__merge_var_with_adjusted(ncIN, 'PRES')

        M_ADJ   = self.__fillnan(ncIN, var + "_ADJUSTED")
        QC      = self.__fillnan(ncIN, var +"_QC")
        QC_adj  = self.__fillnan(ncIN, var +"_ADJUSTED_QC")

        M       = self.__fillnan(ncIN, var )

        Profile     =     M[iProf,:]
        Profile_adj = M_ADJ[iProf,:]
        Pres        =  PRES[iProf,:]
        Qc          =    QC[iProf,:]
        Qc_adj      =QC_adj[iProf,:]
        ncIN.close()

        return Pres, Profile,Profile_adj, Qc.astype(np.int32), Qc_adj.astype(np.int32)

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
        rawPres, rawProfile, rawProfile_adj, rawQc, rawQc_adj = self.read_very_raw(var)

        if read_adjusted:
            rawProfile = rawProfile_adj
            rawQc      = rawQc_adj


        # Elimination of negative pressures or nans
        nanPres = np.isnan(rawPres)
        rawPres[nanPres] = 1 # just for not complaining
        badPres    = (rawPres<=0) | nanPres
        badProfile = np.isnan(rawProfile)
        bad = badPres | badProfile


        Pres    =    rawPres[~bad]
        Profile = rawProfile[~bad]
        Qc      =      rawQc[~bad]

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

    def read(self, var, read_adjusted=True):
        '''

        Reads profile data from file, applies a rarefaction and optionally a filter to the data

        Takes var as string
              read_adjusted as logical
        Returns 3 numpy arrays: Pres, Profile, Qc
        '''
        pres, prof, qc = self.read_raw(var,read_adjusted)
        if pres.size ==0:
            return pres, prof, qc

        if var in ['TEMP','PSAL']:
            good = np.ones_like(pres, dtype=bool)
        else:
            good = (qc==1) | (qc ==2 ) | (qc==5) | (qc==8)
        pres = pres[good]
        prof = prof[good]
        qc   =   qc[good]
        if pres.size ==0:
            return pres, prof, qc

        ii = self.rarefy(pres, 2.0)
        pres = pres[ii]
        prof = prof[ii]
        qc   =   qc[ii]

        if (var=='NITRATE') :
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

        return pres, prof, qc

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
        fig, ax = f.plot(pres,profile,fig,ax, linestyle = 'None',   marker = '.',  color='green')

        '''
        if (fig is None) or (ax is None):
            fig , ax = pl.subplots()
        ax.plot(profile,Pres, **kwargs)
        if not ax.yaxis_inverted(): ax.invert_yaxis()
        ax.grid()
        return fig,ax

    def profiles(self, var):
        return [BioFloatProfile(var, self.time, self.lon, self.lat, self)]

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
            parameterdatamode= INDEX_FILE['parameter_data_mode'][iFile]
            float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')

            if ONLINE_REPO + CORIOLIS_DIR + thefilename == filename :
                return BioFloat(lon,lat,float_time,filename,available_params,parameterdatamode)
        return None

def profile_gen(lon,lat,float_time,filename,available_params,parameterdatamode):

    filename = Path(ONLINE_REPO + CORIOLIS_DIR + filename)
    thefloat = BioFloat(lon,lat,float_time,filename,available_params,parameterdatamode)
    return BioFloatProfile(float_time,lon,lat, thefloat,available_params)

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


    nFiles=INDEX_FILE.size
    selected = []
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        parameterdatamode= INDEX_FILE['parameter_data_mode'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')

        if var is None :
            VarCondition = True
        else:
            VarCondition = var in available_params

        if VarCondition:
            if T.contains(float_time) and region.is_inside(lon, lat):
                selected.append(profile_gen(lon, lat, float_time, filename, available_params,parameterdatamode))

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
    import sys

    f='/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C/CORIOLIS/6902902/SD6902902_456.nc'
    F = BioFloat.from_file(f)

    Pres, values, valuesa, Qc, Qca =F.read_very_raw('TEMP')
    fig,ax =F.plot(Pres,values,'b')

    Pres, values, Qc=F.read('TEMP', read_adjusted=False)
    #ax.plot(values,Pres,'r.')
    #fig.savefig('prova.png')

    var = 'CHLA'
    TI = TimeInterval('20120101','20250225','%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILE_LIST=FloatSelector(var, TI, R)

    nP = len(PROFILE_LIST)
    IND=np.zeros((nP,))
    MIN=np.zeros((nP,))
    for ip, p in enumerate(PROFILE_LIST):
        Pres, Value, Qc= p.read(var, read_adjusted=True)
        if len(Pres)>5:
            MIN[ip]=Pres.min()
            IND[ip]=ip
    sys.exit()

    sum=0
    PROFILE_LIST=FloatSelector(var, TI, R)
    for ip, p in enumerate(PROFILE_LIST):
        if p._my_float.status_var('NITRATE') =='D':
            Pres, Value, Qc= p.read(var, read_adjusted=True)
            print(Pres.max())

        continue
        F=p._my_float
        nP=len(Pres)
        if nP<5 :
            print( "few values for " + F.filename)
            continue
        if Pres[-1]<100:
            print("depth < 100 for "+ F.filename)
            continue

        if Pres.max()< 600:
            print("superficiale", p.ID())
           #sys.exit()

        if F.status_var('NITRATE')=='D':
            #if F.status_var('DOXY') =='R':
            sum+=1
            #print ip, F.filename,'R'

#    filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V5C/FLOAT_BIO/6900807/MR6900807_225.nc"
#    F=BioFloat.from_file(filename)
#    F2=BioFloat.from_file("/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V5C/FLOAT_BIO/6901765/MR6901765_001.nc")




    for p in PROFILE_LIST[:1]:
        PN,N, Qc = p.read(var,read_adjusted=True)
        TheFloat = p._my_float
        PN,N,Qc = TheFloat.read(var,    read_adjusted=True)
        #PS,S,Qc = TheFloat.read('PSAL', read_adjusted=True)
        #PT,T,Qc = TheFloat.read('TEMP', read_adjusted=True)

    wmo_list= get_wmo_list(PROFILE_LIST)
    for wmo in wmo_list:
        sublist = filter_by_wmo(PROFILE_LIST, wmo)
