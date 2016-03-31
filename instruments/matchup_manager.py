import os
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Rectangle

import instruments
import scipy.io.netcdf as NC
import numpy as np

import matchup.matchup
import datetime
import pylab as pl
import matplotlib.patches as mpatches
import seawater as sw
import libxmp
import  libxmp.utils
from libxmp import XMPFiles, consts

class Matchup_Manager():

    '''
    Main class for Float Matchup generation.
    '''
    def __init__(self,timeinterval,INPUTDIR,Outpudir):
        '''
                Outpudir is intended as the outputdir of aveScan,
        point profiles will be produced in outputdir/PROFILES.
        '''
        self.profilingDir=Outpudir
        self.AVE_INPUT_DIR = INPUTDIR
        if os.path.exists(INPUTDIR):
            self.TL = TimeList.fromfilenames(timeinterval, self.AVE_INPUT_DIR,"ave*.nc",'postproc/IOnames.xml')
            self.TI = timeinterval
            All_Med = Rectangle(-6,36,30,46)
            self.PROFILE_LIST = instruments.Selector(None, self.TI, All_Med)
            datetimelist = [f.time for f in self.PROFILE_LIST]
            self.Coupled_List = self.TL.couple_with(datetimelist)
        else:
            print INPUTDIR + " not existing"

    def _dump_punti_for_aveScan(self,Profilelist,filename):
        LINES=[]
        LINES.append('NOME    Longitutine E    Latitudine N \n')
        for p in Profilelist:
            line = "%s\t%g\t%g \n" %(p.name(), p.lon, p.lat)
            LINES.append(line)

        F = open(filename, "w")
        F.writelines(LINES)
        F.close()

    def writefiles_for_profiling(self, filename):
        '''
        Preparation of launch of aveScan.py, in order to generate profiles.

        The file produced - first argument - is a wrapper of aveScan, to call it over times.
        For every launch of aveScan a different punti*.dat files will be used, depending on Biofloats present
        at that time.


        '''

        TMPSDIR =           self.profilingDir + "TMP/"
        PUNTI_DIR =         self.profilingDir + "PUNTI/"

        os.system("mkdir -p "+ PUNTI_DIR)
        JOB_LINES=[]

        JOB_LINES.append("export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc \n")
        JOB_LINES.append("export SUBMASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/POSTPROC/submask.nc \n")
        JOB_LINES.append("cd postproc \n")
        for t in self.Coupled_List:
            Model_time        = t[0]
            INTERESTED_PROFILES = list(set([self.PROFILE_LIST[k] for k in t[1]])) #t[1]


            outpuntifile= PUNTI_DIR + "punti_" + Model_time.strftime("%Y%m%d") + ".dat" #punti_20150416.dat
            self._dump_punti_for_aveScan(INTERESTED_PROFILES, outpuntifile)
            line = 'python aveScan.py '   + \
                ' -l '  + Model_time.strftime("ave.%Y%m%d*")  + \
                ' -i '  + self.AVE_INPUT_DIR  +  \
                ' -t '  + TMPSDIR  + \
                ' -o '  + self.profilingDir  + \
                ' -d VarDescriptorB.xml ' + \
                ' -p ' + outpuntifile + '\n'
            JOB_LINES.append(line)

        F=file(filename,'w')
        F.writelines(JOB_LINES)
        F.close()



    def dumpModelProfiles(self,profilername):
        '''
        The text file provided as argument
         - becomes executable
         - it is launched on cluster
         A series of aveScan.py are executed, once and for all the following analysis.
        '''

        MODEL_PROFILES_DIR =self.profilingDir + "PROFILES/"
        TMPSDIR =           self.profilingDir + "TMP/"
        os.system("mkdir -p " + MODEL_PROFILES_DIR)
        os.system("mkdir -p " + TMPSDIR)
        os.chmod(profilername, 0755)
        os.system(profilername)

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

    def modeltime(self,profile):
        for Model_time,INTERESTED_Indices in self.Coupled_List:
            if profile in [self.PROFILE_LIST[k] for k in INTERESTED_Indices]:
                break
        return Model_time

    def reference_var(self,p,var):
        '''
        returns the reference varname, for a given profile object and
        a model varname
        For BioFloats reference_var(p,'O2o') returns 'DOXY'
        '''
        if isinstance(p, instruments.bio_float.BioFloatProfile):
            return instruments.FLOATVARS[var]
        if isinstance(p, instruments.mooring.MooringProfile):
            return instruments.MOORINGVARS[var]


    def getMatchups(self,Profilelist,nav_lev,model_varname, read_adjusted=True):
        '''
        Float list is a list of Bio_Float objects
        It depends on a user selection in space and time

        SingleFloatMatchup is created with
          - good model values (no nans)
          - good float values (no nans)

        The matchups are couples of values (Model,Ref)
        obtained to a given selection in space and time,
        to be used in statistics.
        Matchups are generated by interpolation from model space to observation space.

        At the moment it can happen that a short model profile is extrapolated over a long float profile.
        Then, a replication of matchups could occur.

        Returns a FloatMatchup istance.
        '''

        Group_Matchup = matchup.matchup.ProfilesMatchup()


        for p in Profilelist:
            Model_time = self.modeltime(p)
            if not self.TI.contains(Model_time) :
                print Model_time.strftime("%Y%m%d-%H:%M:%S is a time not included by profiler")
                continue

            Modelfile = self.profilingDir + "PROFILES/" + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")
            ModelProfile = self.readModelProfile(Modelfile, model_varname, p.name())
            seaPoints = ~np.isnan(ModelProfile)

            if np.isnan(ModelProfile).all() : # potrebbe essere fuori dalla tmask
                print "No model data for (lon,lat) = (%g, %g) " %(p.lon, p.lat)
                continue


            ref_varname = self.reference_var(p, model_varname)
            Pres, Profile, Qc = p.read(ref_varname,read_adjusted)

            MODEL_ON_SPACE_OBS=np.interp(Pres,nav_lev[seaPoints],ModelProfile[seaPoints]).astype(np.float32)

            Matchup = matchup.matchup.ProfileMatchup(MODEL_ON_SPACE_OBS, Profile, Pres, Qc, p)

            Group_Matchup.extend(Matchup)


        return Group_Matchup


    def getFloatMatchups(self,Profilelist,nav_lev,outdir="./"):
        '''
        Float list is a list of Bio_Float objects
        It depends on a user selection in space and time

        SingleFloatMatchup is created with
          - good model values (no nans)
          - good float values (no nans)

        The matchups are couples of values (Model,Ref)
        obtained to a given selection in space and time,
        to be used in statistics.
        Matchups are generated by interpolation from model space to observation space.

        At the moment it can happen that a short model profile is extrapolated over a long float profile.
        Then, a replication of matchups could occur.

        Returns nothing
        '''
        from validation.online.profileplotter import figure_generator
        zlevels_out=np.arange(0,401,5)

        for p in Profilelist:
            Model_time = self.modeltime(p)

            if not self.TI.contains(Model_time) :
                print Model_time.strftime("%Y%m%d-%H:%M:%S is a time not included by profiler")
                continue
            VARLIST = p._my_float.available_params.strip().rsplit(" ")
            VARLIST.remove('PRES')
            Modelfile = self.profilingDir + "PROFILES/" + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")
            read_adjusted = [True,False,False,False,False]

            #density calculator on zlevels_out
            model_varname = 'votemper'
            ref_varname = self.reference_var(p, model_varname)
            ModelProfile = self.readModelProfile(Modelfile, model_varname, p.name())
            seaPoints = ~np.isnan(ModelProfile)
            if np.isnan(ModelProfile).all() : # potrebbe essere fuori dalla tmask
                print "No model data for (lon,lat) = (%g, %g) " %(p.lon, p.lat)
                continue
            Pres, temp, Qc = p.read(ref_varname,read_adjusted[3])
            Temp_out = np.interp(zlevels_out,Pres,temp).astype(np.float32)

            model_varname = 'vosaline'
            ref_varname = self.reference_var(p, model_varname)
            Pres, sal, Qc = p.read(ref_varname,read_adjusted[4])
            sal_out = np.interp(zlevels_out,sal,temp).astype(np.float32)
            density = sw.dens(sal_out,Temp_out,zlevels_out)
            #end density calculator

            correction = [1,1000./density,1,1,1]
            mapgraph = [3,4,5,1,2]
            filename = outdir+"/"+Model_time.strftime('%Y%m%d') +"_"+p.name()
            filepng = filename+".png"
            filenc = filename+".nc"

            f = NC.netcdf_file(filenc, 'w')
            f.createDimension('levels', len(zlevels_out))
            m_array = f.createVariable('lev', 'f', ('levels',))
            m_array[:] = zlevels_out[:]

            fig, axs = figure_generator(p)

            for i,model_varname in enumerate(['P_i','O2o','N3n','votemper','vosaline']):
                ref_varname = self.reference_var(p, model_varname)
                if ref_varname not in VARLIST: continue
                ModelProfile = self.readModelProfile(Modelfile, model_varname, p.name())
                seaPoints = ~np.isnan(ModelProfile)
                Pres, Profile, Qc = p.read(ref_varname,read_adjusted[i])
                if len(Pres) == 0:
                    Pres, Profile, Qc = p.read(ref_varname,not read_adjusted[i])

                print model_varname, len(Profile)
                model_on_common_grid=np.interp(zlevels_out,nav_lev[seaPoints],ModelProfile[seaPoints]).astype(np.float32)
                float_on_common_grid=np.interp(zlevels_out,Pres,Profile).astype(np.float32)
                float_on_common_grid = float_on_common_grid*correction[i]

                Matchup = matchup.matchup.ProfileMatchup(model_on_common_grid, float_on_common_grid, zlevels_out, Qc, p)

                #write on NC file
                name_var = model_varname+"_model"
                m_array = f.createVariable(name_var, 'f', ('levels',))
                m_array[:] = model_on_common_grid[:]
                name_var = model_varname+"_float"
                f_array = f.createVariable(name_var, 'f', ('levels',))
                f_array[:] = float_on_common_grid[:]
                #end write netcdf

                #get subplot
                ax=axs[mapgraph[i]]
                fig, ax = Matchup.plot_subplot(model_varname, fig, ax)
            b_patch = mpatches.Patch(color='red', label='Model')
            g_patch = mpatches.Patch(color='blue', label='Float')
            pl.legend(handles=[b_patch,g_patch], loc=5, borderaxespad=0.)
            pl.savefig(filepng)
            pl.close()

            f.close()
            #img.close_file

            #add metadata to image
            xmpfile = XMPFiles( file_path=filepng, open_forupdate=True )
            xmp = xmpfile.get_xmp()
            if xmp is None:
                xmp = libxmp.XMPMeta()
            xmp.set_property(consts.XMP_NS_DC, 'float', p.name() )
            xmp.set_property(consts.XMP_NS_DC, 'date', Model_time.strftime('%Y%m%d') )
            xmp.set_property(consts.XMP_NS_DC, 'hour', Model_time.strftime('%H:%M:%S') )
            xmp.set_property(consts.XMP_NS_DC, 'position.lat',str(p.lat)+"N")
            xmp.set_property(consts.XMP_NS_DC, 'position.lon',str(p.lon)+"E")
            xmpfile.put_xmp(xmp)
            xmpfile.close_file()


        return
