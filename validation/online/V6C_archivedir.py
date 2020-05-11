import argparse
def argument():
    parser = argparse.ArgumentParser('''
    Prints V6C archive directory for a specific date
    The corresponding ave files in that directory are the best choice for that date,
    and correspond to DU content.
    
    It prints something like that
    analysis/20200505
    forecast/20200507
    ''')
    parser.add_argument(   '--date',"-d",
                                type = str,
                                required = True,
                                help = '20120101')
    
    return parser.parse_args()

args = argument()
from commons import V6C_timing 
print V6C_timing.find_best_datestr(args.date)