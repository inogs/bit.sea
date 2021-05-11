import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Extracts files from archive,
    for V3C version where bio and phys files are stored in the same directory.
    '''
    ,formatter_class=argparse.RawTextHelpFormatter)
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

import datetime
from commons import timerequestors
from dateutil.relativedelta import relativedelta
from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from commons.utils import addsep

now=datetime.datetime.now()
dt=relativedelta(months=1)
datetime_month=now -dt
req=timerequestors.Monthly_req(datetime_month.year, datetime_month.month)


LOC = addsep(args.outdir)
archive_dir= args.arcdir

TI=req.time_interval


for var in ['P_l','O2o','N3n','P_c','N1p','ppn','pH','O3c','CO2airflux', 'pCO2','N4n','O3h','N5s','Z_c']:
    T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/ARCHIVE/',glob_pattern="ave*" +var + ".nc.gz")
    T_bio.extract_analysis(LOC, command="ln -fs $INFILE $OUTFILE", remove_ext=False)



