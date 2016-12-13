import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Extracts statistcs from STAT_PROFILES
    and produces time series plots for run comparison

    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required =True,
                                help = ''' Output dir'''
                                )

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                required =True,
                                help = ''' Input dir'''
                                )

    parser.add_argument(   '--pathtxt', '-t',
                                type = str,
                                required =True,
                                help = ''' Text of path of input files'''
                                )

    parser.add_argument(   '--run', '-r',
                                type = str,
                                required =True,
                                help = ''' Name of run directory'''
                                )

    parser.add_argument(   '--var', '-v',
                                type = str,
                                required =True,
                                help = ''' Name of variable (P_i, ppn, etc.)'''
                                )

    return parser.parse_args()


args = argument()

import os
import os.path as path
import numpy as np
import netCDF4
import pickle

from glob import glob
from datetime import datetime

from commons.mask import Mask
from commons.utils import get_date_string



class CoastEnum:
    coast, open_sea, everywhere = range(3)

    @staticmethod
    def valid(val):
        return val in range(3)

class SubBasinEnum:
    alb, sww, swe, nwm, tyr, adn, ads, aeg, ion, lev, med = range(11)

    @staticmethod
    def valid(val):
        return val in range(11)

class StatEnum:
    mean, std, p25, p50, p75 = range(5)

    @staticmethod
    def valid(val):
        return val in range(5)


SUB_LIST = ['alb', 'sww', 'swe', 'nwm', 'tyr', 'adn', 'ads', 'aeg', 'ion', 'lev', 'med']
COAST_LIST = ['coast', 'open_sea', 'everywhere']
DEPTH_LIST = [0, 50, 150]

maskfile = os.getenv('MASKFILE')
TheMask = Mask(maskfile)

def extract_from_runs(run, varname, subbasin, coast=CoastEnum.open_sea, stat=StatEnum.mean):
    """
    Extracts a time series based on a list of file paths.

    Args:
        - *run*: name of the run
        - *varname*: name of the variable to extract
        - *subbasin*: an element from SubBasinEnum.
        - *coast* (optional): an element from CoastEnum (default:
          CoastEnum.open_sea).
        - *stat* (optional): an element from StatEnum (default: StatEnum.mean).

    Returns: statistics and dates
    """

    MYLIST = list()
    Data_dict = {}
    for dep in DEPTH_LIST: Data_dict[dep] = list()

    run_path = args.indir + run + args.pathtxt
    print(run_path)
    file_list = sorted(glob(run_path + '/*nc'))
    label_list = list()
    for f in file_list:
        #Get date string from file name
        _, ds = get_date_string(path.basename(f))
        #Create datetime object from date string
        dt = datetime.strptime(ds,'%Y%m%d')
        #Append the date to label_list
        label_list.append(dt)
        #Open it with netCDF4
        dset = netCDF4.Dataset(f)
        #Append the variable value to data_list
        for idep,dep in enumerate(DEPTH_LIST):
            depth_index = TheMask.getDepthIndex(dep)
            value = dset[varname][subbasin, coast, depth_index, stat].copy()
            Data_dict[dep].append(value)

        #Close the file
        dset.close()
    print(value)

    LIST = [i for i in range(len(DEPTH_LIST)+1)]
    LIST[0] = label_list
    for idep,dep in enumerate(DEPTH_LIST): LIST[idep+1] = Data_dict[dep]

    return LIST



print(args.var)
for isub in range(len(SUB_LIST)):
    print(SUB_LIST[isub])
    LISTopen = extract_from_runs(args.run, args.var, isub, \
                      coast = CoastEnum.open_sea)
    outfile = args.outdir + args.var + '_' + \
              args.run + '_' + \
              SUB_LIST[isub] + 'open.pkl'
    fid = open(outfile,'wb')
    pickle.dump(LISTopen,fid)
    fid.close()

    LISTcoast = extract_from_runs(args.run, args.var, isub, \
               coast = CoastEnum.coast)
    outfile = args.outdir + args.var + '_' + \
              args.run + '_' + \
              SUB_LIST[isub] + 'coast.pkl'
    fid = open(outfile,'wb')
    pickle.dump(LISTcoast,fid)
    fid.close()

