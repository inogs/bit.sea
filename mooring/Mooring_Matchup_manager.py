import os
from postproc.Timelist import TimeList
from basins.region import Rectangle
from Float.Float_Manager import Time_Interval
import Mooring_Manager
import scipy.io.netcdf as NC
import numpy as np

class Mooring_Matchup_Manager():
    def __init__(self,DATESTART,DATE__END,INPUTDIR,Outpudir):
        '''
        Outpudir is intended as the outputdir of aveScan, 
        point profiles will be produced in outputdir/PROFILES.
        '''
        self.profilingDir=Outpudir
        self.AVE_INPUT_DIR = INPUTDIR
        if os.path.exists(INPUTDIR):
            self.TL = TimeList(DATESTART, DATE__END, self.AVE_INPUT_DIR,"ave*.nc",'../postproc/IOnames.xml')
            self.TI = Time_Interval(DATESTART,DATE__END,'%Y%m%d-%H:%M:%S')
            All_Med = Rectangle(-6,36,30,46)
            self.MOORING_LIST=Mooring_Manager.MooringSelector(None, self.TI, All_Med)
            
            self.profileList=[]
            for m in self.MOORING_LIST:
                self.profileList.extend(m.ProfileID_List())
            datetimelist = [p.time for p in self.profileList]
            self.Coupled_List = self.TL.couple_with(datetimelist)
        else:
            print INPUTDIR + " not existing"
            
    def _dump_punti_for_aveScan(self,Floatlist,filename):
        raise NotImplementedError
    def writefiles_for_profiling(self, filename):
        raise NotImplementedError
    
    def readModelProfile(self, filename,var, wmo):
        '''
        Reads profiles produced by aveScan.py.
        In these files each variable has dimensions (jpk, nProfiles)
        And each profile is identified by the corresponding wmo
        '''
        ncIN = NC.netcdf_file(filename,'r')
        
        M = ncIN.variables[var].data.copy()
        iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
        ncIN.close()
        Profile = M[iProfile,:]
        
        return Profile
    def getModelTime(self,profileid):
        for Model_time,INTERESTED_Indices in self.Coupled_List:
            if profileid in [self.profileList[k] for k in INTERESTED_Indices]:
                break
        return Model_time        
        
    def getMatchups(self,MooringList,nav_lev,model_varname,ref_varname):
        
        for m in MooringList:
            profilelist=m.getProfileList()
            for p in profilelist:
                Model_time = self.getModelTime(p)
                Modelfile = self.profilingDir + "PROFILES/" + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")
                ModelProfile = self.readModelProfile(Modelfile, model_varname, m.name)
                seaPoints = ~np.isnan(ModelProfile)
                
                mooring_profile = m.read_raw(ref_varname,p)