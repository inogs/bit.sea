from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval

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


    parser.add_argument(   '--arcdir', '-adir',
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
end__time=args.endtime
LOC = args.outdir
archive_dir= args.arcdir
# starttime='20160301'
# end__time='20160308'
# LOC = "/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/"
#archive_dir="/pico/home/usera07ogs/a07ogs00/OPA/V4/archive/"
TI=TimeInterval(starttime,end__time,'%Y%m%d')


T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/',glob_pattern="ave*gz")
T_phys= TimeSeries(TI, archive_dir,postfix_dir='OPAOPER_A/'          ,glob_pattern="*gz"   )

T_bio.extract_analysis( LOC + 'output_bio/')
T_phys.extract_analysis(LOC + 'output_phys/');
