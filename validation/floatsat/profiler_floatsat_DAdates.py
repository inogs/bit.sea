#!/usr/bin/env python
# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
# Script to generate profiles of model files in
# the same time and locations where instruments
# such as bioFloats, mooring or vessels have been found.

# When imported, this scripts only defines settings for matchup generation.
from instruments.lovbio_float import FloatSelector

from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Rectangle
import datetime
from commons import genUserDateList as DL

# location of input big ave files, usually the TMP directory.
# ave files are supposed to have N3n, O2o and chl

REFDIR = 'DA_FLOAT_SAT/Winter/'
RUN = 'RUN_SAT_FLOAT_01'

INPUTDIR = '/pico/scratch/userexternal/ateruzzi/' + \
           REFDIR + RUN + '/wrkdir/MODEL/AVE_FREQ_1/'

# output directory, where aveScan.py will be run.

BASEDIR = '/pico/scratch/userexternal/ateruzzi/bit.sea/validation/floatsat/' + \
          REFDIR + RUN + '/PROFILATORE_DAtimes/'


#DATESTART = '20150801'
#DATE__END = '20150907'
DATESTART = '20150101'
DATE__END = '20150203'
DATE__END = '20150131'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
TL_daily = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc",filtervar="N1p")
# TL = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc",filtervar="N1p")
M = TL_daily.getWeeklyList(2)
findfirst = True
for iim in range(len(M)):
    datem = datetime.datetime(M[iim].year,M[iim].month,M[iim].day)
    if (findfirst) & (T_INT.contains(datem)):
        DATESTART_DA = M[iim].string # First Tuesday (first assimilation date)
        findfirst = False

texttime = '-12:00:00'
DAtimelist = DL.getTimeList(DATESTART_DA + texttime,DATE__END + texttime,'days=7')
TL_DA = TimeList(DAtimelist)
TL_DA.inputdir = TL_daily.inputdir
TL_DA.prefix = TL_daily.prefix
TL_DA.filtervar = TL_daily.filtervar

V = TL_daily.getWeeklyList(5)
lastfriday = V[-1]
S = TL_daily.getWeeklyList(6)
firstsaturday = S[0]
T_INTfriday_end = TimeInterval(firstsaturday.string,lastfriday.string,'%Y%m%d')
ALL_PROFILES = FloatSelector(None,T_INTfriday_end, Rectangle(-6,36,30,46))


vardescriptorfile="VarDescriptorB.xml"

#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL_DA,BASEDIR)


    profilerscript = BASEDIR + 'jobProfiler.sh'
    aggregatedir="/pico/scratch/userexternal/ateruzzi/" + \
                 REFDIR + RUN + "/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript, \
                     aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
