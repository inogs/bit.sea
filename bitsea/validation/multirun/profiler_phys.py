#!/usr/bin/env python
# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
# Script to generate profiles of model files in
# the same time and locations where instruments
# such as bioFloats, mooring or vessels have been found.

# When imported, this scripts only defines settings for matchup generation.
#from instruments import instruments
from instruments import lovbio_float
from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
import os
# location of input big ave files, usually the TMP directory.
# ave files are supposed to have N3n, O2o and chl
run = "HC_2017_assw"
run = "DA_Float/RUN_SAT_FLOAT_chl_n/"

INPUTDIR="/gpfs/scratch/userexternal/ateruzzi/" + run + \
    "/wrkdir/MODEL/FORCINGS/"
aggregatedir=INPUTDIR
# output directory, where aveScan.py will be run.
BASEDIR='/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/VALID_float/' + \
    run + '/PROFILATORE_PHYS/'

DATESTART = '20150101'
DATE__END = '20160101'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"T*.nc", prefix="T", hour=0)

import basins.OGS as OGS
ALL_PROFILES = lovbio_float.FloatSelector(None, T_INT, OGS.med)#instruments.getAllProfiles(T_INT)

vardescriptorfile="/galileo/home/userexternal/ateruzzi/bit.sea/validation/multirun/VarDescriptor_valid_online.xml"
#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


    profilerscript = BASEDIR + 'jobProfiler.sh'

    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
