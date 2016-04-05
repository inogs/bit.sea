import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = True,
                                help = 'start date')
    parser.add_argument(   '--endtime','-et',
                                type = str,
                                required = True,
                                help = 'end date')

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')

    parser.add_argument(   '--basedir', '-b',
                                type = str,
                                default = None,
                                required = True,
                                help = '''output directory, where aveScan.py will run.
                                           Usually called PROFILATORE''')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")

    return parser.parse_args()

args = argument()

from commons.time_interval import TimeInterval
from instruments.matchup_manager import Matchup_Manager
import basins.OGS as OGS
from instruments import bio_float
from commons.mask import Mask
from commons.utils import addsep

starttime=args.starttime
end__time=args.endtime
INPUTDIR=args.inputdir
BASEDIR=addsep(args.basedir)


TI=TimeInterval(starttime,end__time,'%Y%m%d')


TheMask=Mask(args.maskfile)

Profilelist=bio_float.FloatSelector(None,TI,OGS.med)
M = Matchup_Manager(TI,INPUTDIR,BASEDIR)

profilerscript = BASEDIR + 'jobProfiler.sh'
M.writefiles_for_profiling(profilerscript) # preparation of data for aveScan
M.dumpModelProfiles(profilerscript) # sequential launch of aveScan

M.getFloatMatchups(Profilelist,TheMask.zlevels,args.outdir)
