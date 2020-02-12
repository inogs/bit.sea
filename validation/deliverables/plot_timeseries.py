# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.2

import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import sys
import numpy as np
import argparse

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
                            required = True,
                            help = 'Input pickle file for coast')

    


    return parser.parse_args()

args = argument()

fid = open(args.open_sea_file)
LIST = pickle.load(fid)
fid.close()
TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_ = LIST


fid = open(args.coast_file)
LIST = pickle.load(fid)
fid.close()
model_coast = LIST[3]


from basins import V2 as OGS
for isub,sub in enumerate(OGS.P):
    print sub.name
    fig, ax = pl.subplots()
    ax.plot(TIMES[:],SAT___MEAN[:,isub],'og',label=' SAT')
    ax.plot(TIMES[:],MODEL_MEAN[:,isub],'-k',label=' RAN')
    ax.plot(TIMES[:],model_coast[:,isub],':k',label=' RAN_coast')
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(0.0, 0.6)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30, fontsize=9)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels,fontsize=10)
    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename=args.outdir+"/"+'chl' + sub.name + ".png"
    pl.savefig(outfilename)
