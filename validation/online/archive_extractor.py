from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from commons.utils import addsep
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''',formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = True,
                                help = 'start date')
    parser.add_argument(   '--endtime','-et',
                                type = str,
                                required = True,
                                help = 'end date')


    parser.add_argument(   '--arcdir', '-a',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")

    return parser.parse_args()

args = argument()

starttime=args.starttime
end__time=args.endtime
LOC = addsep(args.outdir)
archive_dir= args.arcdir

TI=TimeInterval(starttime,end__time,'%Y%m%d')


T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/',glob_pattern="ave*gz")
T_phys= TimeSeries(TI, archive_dir,postfix_dir='OPAOPER_A/'          ,glob_pattern="*gz"   )

T_bio.extract_analysis( LOC + 'output_bio/')
T_phys.extract_analysis(LOC + 'output_phys/');
