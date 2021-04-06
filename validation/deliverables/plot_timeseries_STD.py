# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.2

import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
from matplotlib.ticker import FormatStrFormatter

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-O',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )
    parser.add_argument(   '--open_sea_file', '-o',
                            type = str,
                            required = True,
                            help = 'Input pickle file for open sea')
    parser.add_argument(   '--coast_file', '-c',
                            type = str,
                            required = False,
                            help = 'Input pickle file for coast')

    


    return parser.parse_args()

args = argument()

fid = open(args.open_sea_file)
LIST = pickle.load(fid)
fid.close()
#TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_ = LIST
#TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD = LIST
TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD,CORR = LIST

#fid = open(args.coast_file)
#LIST = pickle.load(fid)
#fid.close()
#model_coast = LIST[3]


from basins import V2 as OGS
for isub,sub in enumerate(OGS.P):
  if (isub != 17):  # DO NOT CONSIDER ATLANTIC SUBBASIN
#   if (sub.name != "adr1"):
    print sub.name
    fig, ax = pl.subplots()
    ax.plot(TIMES,SAT___MEAN[:,isub],'og',label=' SAT')
    ax.fill_between(TIMES,SAT___MEAN[:,isub]-SAT____STD[:,isub],SAT___MEAN[:,isub]+SAT____STD[:,isub],color='palegreen')
    ax.plot(TIMES,MODEL_MEAN[:,isub],'-k',label=' RAN')
    ax.plot(TIMES,MODEL_MEAN[:,isub]-MODEL__STD[:,isub],':k') #,label=' RAN')
    ax.plot(TIMES,MODEL_MEAN[:,isub]+MODEL__STD[:,isub],':k') #,label=' RAN')
#    ax.plot(TIMES,model_coast[:,isub],':k',label=' RAN_coast')
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    ax.tick_params(axis='both', labelsize=12)
    pl.ylim(0.0, np.max(MODEL_MEAN[:,isub]+MODEL__STD[:,isub]) * 1.2 )
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30)
    outfilename=args.outdir+"/"+'chl' + sub.name + "_STD.png"
    pl.savefig(outfilename)
