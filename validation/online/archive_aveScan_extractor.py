import argparse
 
def argument():
    parser = argparse.ArgumentParser(description = '''
    Extracts files from archive
    '''
    ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--date','-d',
                                type = str,
                                required = True,
                                help = 'end date in yyyymmdd format')
    parser.add_argument(   '--arcdir', '-a',
                                type = str,
                                required = True,
                                help = '''Chain archive directory, e.g. /pico/home/usera07ogs/a07ogs00/OPA/V4/archive''')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "Base output directory; inside it output/ will be created.")
 
    return parser.parse_args()
 
args = argument()

from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from commons.utils import addsep

end__time=args.date
OUTDIR = addsep(args.outdir)
archive_dir= "/gpfs/work/OGS_prod_0/OPA/V5C/prod/archive/" #args.arcdir

TI=TimeInterval("20190305",end__time,'%Y%m%d')


T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/STATISTICS/STAT_PROFILES/',glob_pattern="ave*nc")
T_bio.extract_analysis(OUTDIR, command="cp $INFILE $OUTFILE", remove_ext=False)
