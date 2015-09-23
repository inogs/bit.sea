import scipy.io.netcdf as NC
import numpy as np
import os
import postproc
from postproc.Timelist import *
from matchup import *
from region import *
import Float_Manager
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
        for Model_time,INTERESTED_FLOATS in Coupled_List:
            if f.filename in [ff.filename for ff in INTERESTED_FLOATS]:
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
TL = TimeList(DATESTART, DATE__END, INPUTDIR,"ave*.nc",'postproc/IOnames.xml')
TI = Float_Manager.Time_Interval(DATESTART,DATE__END,'%Y%m%d-%H:%M:%S')
R = Rectangle(-6,36,30,46)
ALL_FLOAT_LIST=Float_Manager.FloatSelector(None, TI, R)
datetimelist = [f.time for f in ALL_FLOAT_LIST]
All_Coupled_List = TL.couple_with(datetimelist)
#   ------------- init ------------------------------


MODELVARS={'DOXY':'O2o', \
           'NITRATE':'N3n'}

R1 = Rectangle(-6,10,30,46)
R2 = Rectangle(10,36,30,46)
Layer_1 = Layer(0,50)
Layer_2 = Layer(50,150)
VARLIST=[ 'DOXY','CHL']
SUBLIST = [R1,R2]
LAYERLIST = [ Layer_1, Layer_2]

TI = Float_Manager.Time_Interval("20150903", "20150917", "%Y%m%d")

nVars = len(VARLIST)
nSub  = len(SUBLIST)
nLay  = len(LAYERLIST)
nStats = 5

STATS = np.zeros((nVars,nSub,nLay,nStats),np.float32)


for ivar, ref_varname  in enumerate(VARLIST):
    model_varname = MODELVARS[ref_varname]
    for isub, sub in enumerate(SUBLIST) :
        FLOAT_LIST=Float_Manager.FloatSelector(ref_varname, TI, sub)
        Matchup = getMatchups(FLOAT_LIST,All_Coupled_List,model_varname,ref_varname)
    
        for ilayer, layer in enumerate(LAYERLIST) :
            Mlayer = Matchup.subset(layer)
            
            STATS[ivar, isub, ilayer, 0] = Mlayer.correlation()
            
         

outfilename = "outfile.nc"
ncOUT = NC.netcdf_file(outfilename,'w')
ncOUT.createDimension('nVars', nVars)
ncOUT.createDimension('nSub', nSub)
ncOUT.createDimension('nLay', nLay)
ncOUT.createDimension('nStats',nStats)

ncvar=ncOUT.createVariable('STATS', 'f', (nVars,nSub,nLay,nStats))
ncvar[:] = STATS


ncOUT.close()



# Provare lo stesso per le fisiche
