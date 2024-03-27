import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
import sys
sys.path.append('/g100_scratch/userexternal/camadio0/PPCON/CODE_loss_attention_max_PPCon_CA/')
import matplotlib.pyplot as plt
#from analysis.utils_analysis import from_day_rad_to_day
from dataset_clustering import FloatDataset
import sys
sys.path.append("/g100_scratch/userexternal/camadio0/PPCON/bit.sea/")
from commons.utils import addsep
from commons.time_interval import TimeInterval

print('ok ')
""" Amadio """
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
                            help = 'Input dir')


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
INDIR  = addsep(args.inputdir)
input_file=args.update_file
OUTDIR = addsep(args.outdir)


# check if the file of input start with the string, I splitted the path to avoid 
# error due to the filepath folders name
if input_file.split("/")[-1].startswith('DIFF_floats'):
    print(" .... Dataset is updated")
elif input_file.split("/")[-1].startswith('Float_Index'):
   print  (" .... Dataset is entirely reworked")
else:
   raise ValueError("please check the input file") 



"""end Amadio"""

def make_ds(training_folder,  INDIR):
    """ prepare the directories and call the make_pandas_df function
    """

    if training_folder == "SUPERFLOAT":
        from make_ds_clustering import make_pandas_df

        if not os.path.exists( INDIR  ):
            os.mkdir( INDIR  )

        if not os.path.exists( OUTDIR):
            os.mkdir( OUTDIR  )
        
        print("making ds...")
        make_pandas_df(input_file ,  OUTDIR  + "ds_sf_clustering.csv", INDIR)
        print("superfloat clustering complete dataset created")

    return

print("Preparing file csv for clustering")
make_ds("SUPERFLOAT",  INDIR)  # amadio Inserted INPUDIR
