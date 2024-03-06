import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
import sys
sys.path.append('/g100_scratch/userexternal/camadio0/PPCON/CODE_loss_attention_max_PPCon/')
import matplotlib.pyplot as plt
from analysis.utils_analysis import from_day_rad_to_day
from utils_clustering import make_ds
from dataset_clustering import FloatDataset

""" Amadio """
import argparse
sys.path.append("/g100_work/OGS23_PRACE_IT/COPERNICUS/bit.sea/")
from commons.utils import addsep
def argument():
    parser = argparse.ArgumentParser(description = '''
    inputdir is generally in pwd: V7C/ds/ 
    ds is the folder that contain the SUPERFLOAT directory (e.g. V7C)
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required = True,
                            default = './',
                            help = 'Input dir')
    return parser.parse_args()

args   = argument()
INDIR  = addsep(args.inputdir)
from commons.utils import addsep

"""end Amadio"""

print("Preparing file csv for clustering")
make_ds("SUPERFLOAT", INDIR)  # amadio Inserted INPUDIR

