from commons.time_interval import TimeInterval
from instruments.matchup_manager import Matchup_Manager
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''')


    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = True,
                                help = 'start date')
    parser.add_argument(   '--endtime','-et',
                                type = str,
                                required = True,
                                help = 'end date')

    parser.add_argument(   '--inputdir','-idir',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')

    parser.add_argument(   '--basedir', '-bdir',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")

    parser.add_argument(   '--outdir', '-odir',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")

    return parser.parse_args()

args = argument()

starttime=args.starttime
#starttime='20160301'

end__time=args.endtime
#end__time='20160308'

INPUTDIR=args.inputdir
#INPUTDIR="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/TMP/"

# output directory, where aveScan.py will be run.
BASEDIR=args.basedir
#BASEDIR='/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/PROFILATORE/'

TI=TimeInterval(starttime,end__time,'%Y%m%d')

import basins.OGS as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
from commons.mask import Mask

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')

Profilelist=bio_float.FloatSelector(None,TI,OGS.med)
M = Matchup_Manager(TI,INPUTDIR,BASEDIR)
# # from profiler
#M.writefiles_for_profiling('./jobProfiler.sh') # preparation of data for aveScan
#M.dumpModelProfiles('./jobProfiler.sh') # sequential launch of aveScan
# end from profiler
M.getFloatMatchups(Profilelist,TheMask.zlevels,args.outdir)
