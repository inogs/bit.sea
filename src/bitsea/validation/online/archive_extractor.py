import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Extracts files from archive,
    for V3C version where bio and phys files are stored in the same directory.
    '''
    ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-et',
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
    parser.add_argument(   '--type', 
                                type = str,
                                choices = ['analysis','forecast'],
                                required = True)

    return parser.parse_args()

args = argument()

from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from commons.utils import addsep

starttime=args.starttime
end__time=args.endtime
LOC = addsep(args.outdir)
archive_dir= args.arcdir

TI=TimeInterval(starttime,end__time,'%Y%m%d')

if args.type=='analysis':
    for var in ['P_l','O2o','N3n','vosaline','votemper','EIR','P_c','pH','POC']:
        T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/ARCHIVE/',glob_pattern="ave*" +var + ".nc.gz")
        T_bio.extract_analysis( LOC + 'output/')
#     T_phys= TimeSeries(TI, archive_dir,postfix_dir='OPAOPER_A/'          ,glob_pattern="*T.nc"   )
#     T_phys.extract_analysis(LOC + 'output_phys_ingv/', command="cp $INFILE $OUTFILE", remove_ext=False);


if args.type =='forecast':

    for var in ['P_l','O2o','N3n','vosaline','votemper','EIR','P_c','pH','POC']:
       T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/ARCHIVE/',glob_pattern="ave*" +var + ".nc.gz")
       T_bio.extract_simulation(LOC + 'output/')
       T_bio.extract_forecast(  LOC + 'output/')

#     T_phys_s= TimeSeries(TI, archive_dir,postfix_dir='OPAOPER_A/'          ,glob_pattern="*T.nc" )
#     T_phys_f= TimeSeries(TI, archive_dir,postfix_dir='OPAOPER_F/'          ,glob_pattern="*T.nc" )
#     T_phys_s.extract_simulation(LOC + 'output_phys_ingv/', command="cp $INFILE $OUTFILE", remove_ext=False);
#     T_phys_f.extract_forecast(  LOC + 'output_phys_ingv/', command="cp $INFILE $OUTFILE", remove_ext=False);
