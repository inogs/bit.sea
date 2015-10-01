import scipy.io.netcdf as NC
import numpy as np
import os
from postproc.Timelist import *
from matchup.matchup import *
from basins.region import *
import Float.Float_Manager as Float_Manager
from shared_data import *


def getModelProfile(filename,var, wmo):
    ncIN = NC.netcdf_file(filename,'r')
    
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile,:]
    
    return Profile

def getMatchups(FloatList,Coupled_List,model_varname,ref_varname):
    ''' 
    Float list is a list of Bio_Float objects 
    It depends on a user selection in space and time
    
    Coupled list is a list of tuples (datetime object, [list of integers] ) 
    as returned by Timelist.couple_with()
    It depends on a selection used to launch the jobProfiler.
    '''

    
    Group_Matchup = FloatMatchup()
    
    
    for f in FloatList:
        for Model_time,INTERESTED_Indices in Coupled_List:
            if f.filename in [ALL_FLOAT_LIST[ii].filename for ii in INTERESTED_Indices]:
                break
        Modelfile = MODEL_PROFILES_DIR + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")
        ModelProfile = getModelProfile(Modelfile, model_varname, f.wmo)
        seaPoints = ~np.isnan(ModelProfile)
                
        if np.isnan(ModelProfile).all() : # potrebbe essere fuori dalla tmask         
            print "No model data for (lon,lat) = (%g, %g) " %(f.lon, f.lat)
            continue
        
        FloatPres, FloatProfile = f.read(ref_varname)
        
        MODEL_ON_SPACE_OBS=np.interp(FloatPres,nav_lev[seaPoints],ModelProfile[seaPoints]).astype(np.float32)
                
        Matchup = SingleFloatMatchup(MODEL_ON_SPACE_OBS, FloatProfile, FloatPres,f)
        
        Group_Matchup.extend(Matchup)

        
    return Group_Matchup  


maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

# getting Coupled list for all the data
# ----------------  init ------------------------------
TL = TimeList(DATESTART, DATE__END, INPUTDIR,"ave*.nc",'../postproc/IOnames.xml')
TI = Float_Manager.Time_Interval(DATESTART,DATE__END,'%Y%m%d-%H:%M:%S')
R = Rectangle(-6,36,30,46)
ALL_FLOAT_LIST=Float_Manager.FloatSelector(None, TI, R)
datetimelist = [f.time for f in ALL_FLOAT_LIST]
All_Coupled_List = TL.couple_with(datetimelist)
#   ------------- init ------------------------------


MODELVARS={'DOXY':'O2o', \
           'NITRATE':'N3n',  \
           'CHLA':'P_i'}
           

R1 = Rectangle(-6,10,30,46)
R2 = Rectangle(10,36,30,46)


VARLIST=[ 'DOXY','CHLA']
SUBLIST = [R1,R2]
import basins.OGS as OGS
SUBLIST = [sub for sub in OGS.P]

LAYERLIST = [Layer(0,10), Layer(10,50), Layer(50,100), Layer(100,150), Layer(150,300),Layer(300,600), Layer(600,1000)]

TI = Float_Manager.Time_Interval("20150903", "20150917", "%Y%m%d")

nVars = len(VARLIST)
nSub  = len(SUBLIST)
nLay  = len(LAYERLIST)

METRICS_NAMES=['number of data values',\
               'mean of product',\
               'mean of reference',\
               'mean squared error',\
               'variance of product',\
               'variance of reference']

nMetrics = len(METRICS_NAMES)

STATS = np.ones((nVars,nSub,nLay,nMetrics),np.float32)*np.nan


for ivar, ref_varname  in enumerate(VARLIST):
    model_varname = MODELVARS[ref_varname]
    print model_varname
    for isub, sub in enumerate(SUBLIST) :
        print sub
        FLOAT_LIST=Float_Manager.FloatSelector(ref_varname, TI, sub)
        Matchup = getMatchups(FLOAT_LIST,All_Coupled_List,model_varname,ref_varname)
    
        for ilayer, layer in enumerate(LAYERLIST) :
            print layer
            Mlayer = Matchup.subset(layer)
            
            n = Mlayer.number()
            if n>0:
                STATS[ivar, isub, ilayer, 0] = n 
                STATS[ivar, isub, ilayer, 1] = Mlayer.Model.mean()
                STATS[ivar, isub, ilayer, 2] = Mlayer.Ref.mean()
                STATS[ivar, isub, ilayer, 3] = Mlayer.MSE()
                STATS[ivar, isub, ilayer, 4] = np.median(Mlayer.Model)
                STATS[ivar, isub, ilayer, 5] = np.median(Mlayer.Ref)
         

outfilename = "outfile.nc"
ncOUT = NC.netcdf_file(outfilename,'w')
ncOUT.createDimension('nVars', nVars)
ncOUT.createDimension('nSub', nSub)
ncOUT.createDimension('nLay', nLay)
ncOUT.createDimension("metrics",nMetrics     )


ncvar=ncOUT.createVariable('METRICS', 'f', ('nVars','nSub','nLay','metrics'))
ncvar[:] = STATS


ncOUT.close()



# Provare lo stesso per le fisiche
