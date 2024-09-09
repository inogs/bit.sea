import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--reportdir', '-rep',
                                type = str,
                                default = None,
                                required = True,
                                help = '''reports directoryes e.g. OUTPUTS''')
    return parser.parse_args()

args = argument()

import pandas as pd 
import numpy as np
import glob
from commons_ import col_to_dt
from commons.utils import addsep

OUTDIR      = addsep(args.reportdir)

#OUTDIR='/g100/home/userexternal/camadio0/float_preproc/TimeSeries_plots/Superfloat_oxy/OUTPUTS/'

filedir = glob.glob( OUTDIR  + "*_report.csv")

li = []

datelist=[]
for filename in filedir:
    df = pd.read_csv(filename, index_col=0, header=0)
    li.append(df)
    Dataitem = filename.replace(OUTDIR,"").replace("_report.csv","")
    datelist.append(Dataitem)

datelist.sort(reverse=False)
frame = pd.concat(li, axis=0, ignore_index=True)

frame.to_csv(OUTDIR +'Report_'+ datelist[0] + '_' + datelist[-1]+ ".csv"  )


