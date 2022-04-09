import argparse
def argument():
    parser = argparse.ArgumentParser('''
    Prints V8C archive info for a specific date
    The corresponding files are the best choice for that date,
    and correspond to DU content.
    
    It prints prefixes of files, something like that
    1. with --phys : forecast/20200511/CMCC_PHYS/mfs_eas5-20200511-20200512-f-              ( + T.nc)
    2. with --bgc  : forecast/20200511/POSTPROC/AVE_FREQ_1/ARCHIVE/ave.20200512-12:00:00.   ( + var.nc.gz)
    3. with --dir  : analysis/20200505
    ''')
    parser.add_argument(   '--date',"-d",
                                type = str,
                                required = True,
                                help = '20120101')
    parser.add_argument(   '--phys',
                                action='store_true',
                                required = False,
                                help = 'enables physical forcing print')
    parser.add_argument(   '--bgc',
                                action='store_true',
                                required = False,
                                help = 'enables biogeochemistry print')
    parser.add_argument(   '--dir',
                                action='store_true',
                                required = False,
                                help = 'enables directory print')
    parser.add_argument(   '--maps',
                                action='store_true',
                                required = False,
                                help = 'enables ordered printout of maps.tar')
    parser.add_argument(   '--analysis_already_run',
                                action='store_true',
                                required = False,
                                help = 'True if we are sure that analysis of today is already archived')
    parser.add_argument(   '--forecast_already_run',
                                action='store_true',
                                required = False,
                                help = 'True if we are sure that forecast of today is already archived')    
    
    return parser.parse_args()

args = argument()
from commons import V8C_timing as timing
if args.phys:
    print(timing.find_best_forcing(args.date,sure_an_already_run=args.analysis_already_run, sure_fc_already_run=args.forecast_already_run))
if args.bgc:
    print(timing.find_best_bgc(args.date))
if args.dir:
    print(timing.find_best_dir(args.date))
if args.maps:
    for thedir in timing.list_for_maps(args.date, sure_an_already_run=args.analysis_already_run, sure_fc_already_run=args.forecast_already_run):
        print(thedir + "/POSTPROC/AVE_FREQ_1/maps.tar")
