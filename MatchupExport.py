# This script is meant to generate matchups for
# all the Med and export them into text files, without any local calculation.
# All setup data are imported from profiler.py
# There are some restrictions: use of only BioFloats, only chl


from instruments.bio_float import *
from  basins.region import Rectangle
from profiler import * 

OUTDIR  = "/pico/scratch/userexternal/gbolzon0/E-AIMS-GP/Matchup/CNTRL/"




M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)

maskfile    = os.getenv("MASKFILE");
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

All_Med = Rectangle(-6,36,30,46)
#a=TimeInterval('20130301','20130501','%Y%m%d')
Profilelist=FloatSelector('CHLA',T_INT,All_Med)

Matchup = M.getMatchups(Profilelist,nav_lev,'P_i')

prefix=""

Matchup.export(OUTDIR,prefix)




