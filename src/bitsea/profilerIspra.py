#!/usr/bin/env python
# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
# Script to generate profiles of model files in
# the same time and locations where instruments
# such as bioFloats, mooring or vessels have been found.

# When imported, this scripts only defines settings for matchup generation.
from static import Ispra_reader

from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Rectangle

# location of input big ave files, usually the TMP directory.
# ave files are supposed to have N3n, O2o and chl
RUN = '16'
INPUTDIR="/pico/scratch/userexternal/ateruzzi/DA_COAST_" + RUN + "/wrkdir/MODEL/AVE_FREQ_1"

# output directory, where aveScan.py will be run.
BASEDIR='/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/PROFILATORE_' + RUN + '/'

DATESTART = '20130101-00:00:00'
DATE__END = '20140101-00:00:00'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"*.nc",filtervar='N1p')

Reg= Rectangle(-6,36,30,47)
N = Ispra_reader.IspraReader()
ISPRA_PROFILES = N.Selector(None, T_INT, Reg)

vardescriptorfile="VarDescriptorB.xml" # in bit.sea/postproc
#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ISPRA_PROFILES,TL,BASEDIR)


    profilerscript = BASEDIR + 'jobProfiler_' + RUN + '.sh'

    AGGDIR="/pico/scratch/userexternal/ateruzzi/DA_COAST_" + RUN + "/wrkdir//POSTPROC/output/AVE_FREQ_1/TMP/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=AGGDIR) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
