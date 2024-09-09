import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates a png and a NetCDF file for each float
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
                                help = 'input dir with bio ave 3D')
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
                                required = False,
                                help = "directory to store png and NetCDF files")
    parser.add_argument(   '--xml', '-x',
                                type = str,
                                default = 'VarDescriptor_valid_online.xml',
                                required = False,
                                help = "Var descriptor file, the default is bit.sea/postproc/VarDescriptor_valid_online.xml")


    return parser.parse_args()

args = argument()
import matplotlib
matplotlib.use('Agg')
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from instruments.matchup_manager import Matchup_Manager
import basins.OGS as OGS
from instruments import superfloat as bio_float
from commons.mask import Mask
from commons.utils import addsep
from datetime import timedelta

starttime=args.starttime
end__time=args.endtime
INPUTDIR=addsep(args.inputdir)
BASEDIR=addsep(args.basedir)


TI=TimeInterval(starttime,end__time,'%Y%m%d')
TI.end_time = TI.end_time + timedelta(1)

TheMask=Mask(args.maskfile)

Profilelist=bio_float.FloatSelector(None,TI,OGS.med)
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc",filtervar="P_l")
M = Matchup_Manager(Profilelist,TL,BASEDIR)

profilerscript = BASEDIR + 'jobProfiler.sh'
M.writefiles_for_profiling(args.xml, profilerscript, aggregatedir=INPUTDIR) # preparation of data for aveScan
M.dumpModelProfiles(profilerscript) # sequential launch of aveScan

if args.outdir:
    M.getFloatMatchups(Profilelist,TheMask.zlevels,args.outdir)
