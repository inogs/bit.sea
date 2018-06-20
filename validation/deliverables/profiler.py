#!/usr/bin/env python
# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
# Script to generate profiles of model files in
# the same time and locations where instruments
# such as bioFloats, mooring or vessels have been found.

# When imported, this scripts only defines settings for matchup generation.
#from instruments import instruments
from instruments import lovbio_float
from static.Massimili_reader import MassimiliReader
from instruments.matchup_manager import Matchup_Manager
from commons.Timelist import TimeInterval, TimeList
import basins.V2 as OGS
# location of input big ave files, usually the TMP directory.
# ave files are supposed to have N3n, O2o and chl
INPUTDIR="/pico/scratch/userexternal/lmariott/FDA_all2015_newstd_3days/wrkdir/MODEL/AVE_FREQ_1/"
aggregatedir="/pico/scratch/userexternal/lmariott/FDA_all2015_newstd_3days/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/"
# output directory, where aveScan.py will be run.
BASEDIR='/pico/scratch/userexternal/gbolzon0/MASSIMILI/FDA_all2015_newstd_3days/PROFILATORE/'

DATESTART = '20150101-00:00:00'
DATE__END = '20160101-00:00:00'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc", filtervar="N1p")


N = MassimiliReader()
ALL_PROFILES = N.Selector(None, T_INT, OGS.med) #lovbio_float.FloatSelector(None, T_INT, OGS.med)#instruments.getAllProfiles(T_INT)

vardescriptorfile="VarDescriptorB.xml"
#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


    profilerscript = BASEDIR + 'jobProfiler.sh'

    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
