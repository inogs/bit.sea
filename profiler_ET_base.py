#!/usr/bin/env python
# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
# Script to generate profiles of model files in
# the same time and locations where instruments
# such as bioFloats, mooring or vessels have been found.

# When imported, this scripts only defines settings for matchup generation.
#from instruments import instruments
from instruments import optbio_float_2019
from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
# location of input big ave files, usually the TMP directory.
# ave files are supposed to have N3n, O2o and chl
INPUTDIR="/gpfs/scratch/userexternal/eterzic0/BIOPTIMOD/<DIR>/<SUBDIR>/INPUT/"#"/gpfs/scratch/userexternal/plazzari/REA_16_T0/TEST_01/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/"

# output directory, where aveScan.py will be run.
BASEDIR="/gpfs/scratch/userexternal/eterzic0/BIOPTIMOD/<DIR>/<SUBDIR>/PROFILATORE/"

DATESTART = '20120101-00:00:00'
DATE__END = '20180101-00:00:00'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc",filtervar="Ed380f")

import basins.OGS as OGS
ALL_PROFILES = optbio_float_2019.FloatSelector(None, T_INT, OGS.med)#instruments.getAllProfiles(T_INT)

vardescriptorfile="/galileo/home/userexternal/eterzic0/BIOPTIMOD/<DIR>/<SUBDIR>/bit.sea/Vardescriptor_bioptimod.xml"
#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


    profilerscript = BASEDIR + 'jobProfiler.sh'
    aggregatedir="/gpfs/scratch/userexternal/eterzic0/BIOPTIMOD/<DIR>/<SUBDIR>/INPUT/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript,aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
