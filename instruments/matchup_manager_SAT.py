import os,sys
from instrument import Profile
import scipy.io.netcdf as NC
import numpy as np
import matchup.matchup
import pylab as pl
import all_instruments
import postproc


class Matchup_Manager():

    '''
    Main class for Float Matchup generation.
    '''
    def __init__(self,PROFILE_LIST,Model_Timelist,Outpudir):
        '''
            Model_Timelist is a Timelist object
                Outpudir is intended as the outputdir of aveScan,
        point profiles will be produced in outputdir/PROFILES.
        '''
        self.profilingDir=Outpudir
        self.AVE_INPUT_DIR = Model_Timelist.inputdir

        self.TL = Model_Timelist#TimeList.fromfilenames(timeinterval, self.AVE_INPUT_DIR,"ave*.nc")
        self.TI = Model_Timelist.timeinterval
        self.PROFILE_LIST = PROFILE_LIST
        datetimelist = [f.time for f in self.PROFILE_LIST]
        self.Coupled_List = self.TL.couple_with(datetimelist)


    def _dump_punti_for_aveScan(self,Profilelist,filename):
        LINES=[]
        LINES.append('NOME    Longitutine E    Latitudine N \n')
        for p in Profilelist:
            line = "%s\t%g\t%g \n" %(p.ID(), p.lon, p.lat)
            LINES.append(line)

        F = open(filename, "w")
        F.writelines(LINES)
        F.close()

    def writefiles_for_profiling(self, vardescriptor,filename, aggregatedir="./", ionamesfile="IOnames.xml"):
        '''
        Preparation of launch of aveScan.py, in order to generate profiles.
        Arguments
        * vardescriptor * a file in postproc/ directory
        * filename     *  the output file, is a wrapper of aveScan, to call it over times.
        * aggregatedir * is the path of the 2nd directory where aveScan will search files

        For every launch of aveScan a different punti*.dat files will be used, depending on Biofloats present
        at that time.


        '''

        TMPSDIR =           self.profilingDir + "TMP/"
        PUNTI_DIR =         self.profilingDir + "PUNTI/"

        os.system("mkdir -p "+ PUNTI_DIR)
        JOB_LINES=[]
        JOB_LINES.append("cd " + os.path.abspath(postproc.__path__[0]) + " \n")
        for t in self.Coupled_List:
            Model_time        = t[0]
            Coupled_indexes   = t[1]
            INTERESTED_PROFILES = []
            for k in Coupled_indexes:
                p=self.PROFILE_LIST[k]
                if not p in INTERESTED_PROFILES:
                    INTERESTED_PROFILES.append(p)
            #INTERESTED_PROFILES = list(set([self.PROFILE_LIST[k] for k in t[1]])) #t[1]

            if self.TL.filtervar is None:
                filtervarline=""
            else:
                filtervarline= ' -f '  + self.TL.filtervar
            outpuntifile= PUNTI_DIR + "punti_" + Model_time.strftime("%Y%m%d") + ".dat" #punti_20150416.dat
            self._dump_punti_for_aveScan(INTERESTED_PROFILES, outpuntifile)
            line = 'python aveScan.py '   + \
                ' -l '  + self.TL.prefix + Model_time.strftime("%Y%m%d*")  + \
                ' -i '  + self.AVE_INPUT_DIR  +  \
                ' -a '  + aggregatedir       +  \
                filtervarline                +  \
                ' -t '  + TMPSDIR  + \
                ' -o '  + self.profilingDir  + \
                ' -d '  + vardescriptor      + \
                ' --ionames ' + ionamesfile  + \
                ' -p ' + outpuntifile + '\n'
            JOB_LINES.append(line)
            JOB_LINES.append("if [ $? -ne 0 ] ; then echo ERROR in aveScan; exit 1  ; fi \n")

        F=file(filename,'w')
        F.writelines(JOB_LINES)
        F.close()



    def dumpModelProfiles(self,profilername, erase=False):
        '''
         A series of aveScan.py are executed, once and for all the following analysis.

        Arguments:
        * profilername * string, text file
        * erase        * logica, if True PROFILES and TMP/ are erased and re-created

        profilername becoomes executable and itit is launched on cluster

        '''
        MODEL_PROFILES_DIR =self.profilingDir + "PROFILES/"
        TMPSDIR =           self.profilingDir + "TMP/"
        if erase:
            os.system("rm -rf " + MODEL_PROFILES_DIR)
            os.system("rm -rf " + TMPSDIR)
            os.system("mkdir " + MODEL_PROFILES_DIR)
            os.system("mkdir " + TMPSDIR)
        os.chmod(profilername, 0755)
        os.system(profilername)
    @staticmethod
    def readModelProfile(filename,var, wmo):
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

    def modeltime(self,profile):
        '''
        Argument:
           * profile * a profile object.
                       E.g. an element of list provided by FloatSelector(), or MooringSelector()
        Returns:
           * datetime object * , the time of model corresponding to the selected profile

        '''
        assert isinstance(profile, Profile)

        for Model_time,INTERESTED_Indices in self.Coupled_List:
            if profile in [self.PROFILE_LIST[k] for k in INTERESTED_Indices]:
                break
        else:
            raise ValueError("No model profile corresponds to %s. Please check time interval in your profiler." % ( profile.ID(), ) )
            return None
        return Model_time

    def reference_var(self,p,var):
        '''
        returns the reference varname, for a given profile object and
        a model varname
        For BioFloats reference_var(p,'O2o') returns 'DOXY'
        '''
        if isinstance(p, all_instruments.optbio_float.BioFloatProfile):
            return all_instruments.FLOAT_OPT_VARS[var]
        if isinstance(p, all_instruments.optbio_float_2019.BioFloatProfile):
            return all_instruments.FLOAT_OPT_BIOPTIMOD[var]
        if isinstance(p, all_instruments.superfloat.BioFloatProfile):
            return all_instruments.FLOATVARS[var]
        if isinstance(p, all_instruments.bio_float.BioFloatProfile):
            return all_instruments.FLOATVARS[var]
        if isinstance(p, all_instruments.lovbio_float.BioFloatProfile):
            return all_instruments.LOVFLOATVARS[var]
        if isinstance(p, all_instruments.mooring.MooringProfile):
            return all_instruments.MOORINGVARS[var]


    def getMatchups(self,Profilelist,nav_lev,model_varname, read_adjusted=True, interpolation_on_Float=True):
        '''
        Arguments
        * Profilelist            * a list of profile objects, generated by a user selection in space and time
        * nav_lev                * array of model depths
        * model_varname          * string
        * read_adjusted          * boolean, if True ADJUSTED variable in NetCDF files will be read
        * interpolation_on_Float * boolean,
                                   if True,  matchups are generated by interpolation from model space to observation space.
                                   if False, matchups are generated by interpolation from observation space to model


        The matchups are couples of values (Model,Ref)
        obtained to a given selection in space and time,
        to be used in statistics.


        At the moment it can happen that a short model profile is extrapolated over a long float profile.
        Then, a replication of matchups could occur.

        Returns a FloatMatchup istance.
        '''

        Group_Matchup = matchup.matchup.ProfilesMatchup()


        for p in Profilelist:
            assert p in self.PROFILE_LIST
            Model_time = self.modeltime(p)
            if not self.TI.contains(Model_time) :
                print Model_time.strftime("%Y%m%d-%H:%M:%S is a time not included by profiler")
                continue

            Modelfile = self.profilingDir + "PROFILES/" + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")
            ModelProfile_aux = self.readModelProfile(Modelfile, model_varname, p.ID())
#BIOPTIMOD            
            ModelProfile = ModelProfile_aux[0]
#BIOPTIMOD            
            seaPoints = ~np.isnan(ModelProfile)

            if np.isnan(ModelProfile).all() : # potrebbe essere fuori dalla tmask
                print "No model data for ", p.ID()
                continue


            ref_varname = self.reference_var(p, model_varname)
            if isinstance(p, (all_instruments.superfloat.BioFloatProfile, all_instruments.optbio_float.BioFloatProfile, all_instruments.optbio_float_2019.BioFloatProfile)):
                Pres, Profile, Qc = p.read(ref_varname)
            else:
                Pres, Profile, Qc = p.read(ref_varname,read_adjusted)

            if interpolation_on_Float:
                MODEL_ON_SPACE_OBS=np.interp(Pres,nav_lev[seaPoints],ModelProfile[seaPoints]).astype(np.float32)

                Matchup = matchup.matchup.ProfileMatchup(MODEL_ON_SPACE_OBS, Profile, Pres, Qc, p)
            else:
                OBS_ON_SPACE_MODEL=np.interp(nav_lev[seaPoints], Pres, Profile)
                QC_ON_SPACE_MODEL = np.interp(nav_lev[seaPoints], Pres, Qc)
                Matchup = matchup.matchup.ProfileMatchup(ModelProfile[seaPoints], OBS_ON_SPACE_MODEL, nav_lev[seaPoints], QC_ON_SPACE_MODEL, p)

            #Group_Matchup.extend(Matchup)

        return Matchup
        #return Group_Matchup


    def getFloatMatchups(self,Profilelist,nav_lev,outdir="./"):
        '''
        Dumps an image png file and a NetCDF file for each profile
        The filenames refers to time and float wmo.
        Both files contain matchups about chl,O2o,N3n,temperature, salinity


        Matchups are provided at fixed depths: every 5 meters from 0 to 400m,
        because in biofloats physical and biological variables do not have the same sampling.

        Arguments:
        * Profilelist * is provided by FloatSelector
        * nav_lev * is the model level
        * outdir * optional is the output location


        The matchups are couples of values (Model,Ref)
        obtained to a given selection in space and time,
        to be used in statistics.


        At the moment it can happen that a short model profile is extrapolated over a long float profile.
        Then, a replication of matchups could occur.

        Returns nothing
        '''
        from validation.online.profileplotter import figure_generator, ncwriter#, add_metadata
        zlevels_out=np.arange(0,401,5)
        MODELVARLIST=['P_l','O2o','N3n','votemper','vosaline','EIR','POC',"P_c", "pH"]
        plotvarname = [r'Chl $[mg/m^3]$',
                       r'Oxy $[mmol/m^3]$',
                       r'Nitr $[mmol/m^3]$',
                       r'Temp $[^\circ C]$',
                       'Sal [psu]',
                       r'PAR $[\mu E/m^2 s]$',
                       r'POC $[mg/m^3]$',
                       'PhytoC $[mg/m^3]$',
                       'pH'
                        ]


        mapgraph = [5,6,7,1,2,8,9,3,4]

        for p in Profilelist:
            Model_time = self.modeltime(p)

            if Model_time is None :
                print p.time.strftime("%Y%m%d-%H:%M:%S is a time not included by profiler")
                continue
            else:
                if not self.TI.contains(Model_time) :
                    print Model_time.strftime("%Y%m%d-%H:%M:%S is a time not included by profiler")
                    continue
            VARLIST = p._my_float.available_params.strip().rsplit(" ")
            if "PRES" in VARLIST: VARLIST.remove('PRES')
            Modelfile = self.profilingDir + "PROFILES/" + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")


            model_varname = 'votemper'
            ref_varname = self.reference_var(p, model_varname)
            ModelProfile = self.readModelProfile(Modelfile, model_varname, p.ID())
            if np.isnan(ModelProfile).all() : # potrebbe essere fuori dalla tmask
                print "No model data for (lon,lat) = (%g, %g) " %(p.lon, p.lat)
                continue

            filename = outdir+"/"+Model_time.strftime('%Y%m%d') +"_"+p.name()
            ncOUT, model_handlers, float_handlers =ncwriter(filename+".nc", zlevels_out,p)

            fig, axs = figure_generator(p)
            #pl.rc('text', usetex=True)

            for i,model_varname in enumerate(MODELVARLIST):
                ref_varname = self.reference_var(p, model_varname)
                if ref_varname not in VARLIST: continue
                ModelProfile = self.readModelProfile(Modelfile, model_varname, p.ID())
                seaPoints = ~np.isnan(ModelProfile)
                Pres, Profile, Qc = p.read(ref_varname)
                if len(Pres) == 0:
                    Pres, Profile, Qc = p.read(ref_varname)

                print model_varname, len(Profile)
                if len(Profile) < 2 : continue
                model_on_common_grid=np.interp(zlevels_out,nav_lev[seaPoints],ModelProfile[seaPoints]).astype(np.float32)
                float_on_common_grid=np.interp(zlevels_out,Pres,Profile).astype(np.float32)
                if model_varname=='POC':
                    float_on_common_grid = float_on_common_grid *  52779.37 - 3.57 # Bellacicco 2019
                if model_varname == 'P_c':
                     bbp470 = float_on_common_grid * ( 470.0/ 700)** 0.78# [m-1]
                     float_on_common_grid = 12128 * bbp470 + 0.59
                if model_varname == "PAR":
                    sec = p.time.hours*3600 + p.time.min*60
                    POSITIVE_VALUE = np.cos(2*np.pi*sec/86400. -np.pi )
                    model_on_common_grid = np.max(0.001, np.pi*POSITIVE_VALUE  *     model_on_common_grid )



                Matchup = matchup.matchup.ProfileMatchup(model_on_common_grid, float_on_common_grid, zlevels_out, Qc, p)

                model_handlers[i][:] = model_on_common_grid[:] #write on NC file
                float_handlers[i][:] = float_on_common_grid[:]


                ax=axs[mapgraph[i]] #get subplot
                fig, ax = Matchup.plot_subplot(plotvarname[i], fig, ax)

            ncOUT.close()
            pngfile = filename + ".png"
            fig.savefig(pngfile)
            pl.close(fig)
            #add_metadata(pngfile, p)


        return
