import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    inputdir is generally in pwd: ONLINE_REPO
    ds is the folder that contain the SUPERFLOAT directory (e.g. V7C)
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required = True,
                            default = './',
                            help = 'an ONLINE REPO dir')


    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required = True,
                            default = './clustering/',
                            help = 'Outdir dir')

    parser.add_argument(   '--update_file','-u',
                                type = str,
                                required = True,
                                default = '',
                                help = ''' or DIFF_floats[DATE].txt files or Float_Index.txt   ''')

    return parser.parse_args()

args   = argument()

import os
import pandas as pd
import torch
from torch.utils.data import DataLoader
from dataset_clustering import FloatDataset
from commons.utils import addsep
from make_ds_clustering import make_pandas_df


INDIR  = addsep(args.inputdir)
input_file=args.update_file
OUTDIR = addsep(args.outdir)


# check if the file of input start with the string, I splitted the path to avoid 
# error due to the filepath folders name
if input_file.split("/")[-1].startswith('DIFF_floats'):
    pass
elif input_file.split("/")[-1].startswith('Float_Index'):
    print  (" .... Dataset is entirely reworked")
else:
    raise ValueError("Input file name must start with 'DIFF_floats' or 'Float_Index' ")


print("Preparing file csv for clustering")
make_pandas_df(input_file ,  OUTDIR  + "ds_sf_clustering.csv", INDIR)
print("superfloat clustering complete dataset created")





