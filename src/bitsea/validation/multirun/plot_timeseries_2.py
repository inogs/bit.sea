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
    parser.add_argument(   '--open_sea_fileV2', '-oV2',
                            type = str,
                            required = True,
                            help = 'Input pickle file for open sea')
    parser.add_argument(   '--coast_fileV2', '-cV2',
                            type = str,
                            required = True,
                            help = 'Input pickle file for coast')
    parser.add_argument(   '--open_sea_fileV3', '-oV3',
                            type = str,
                            required = True,
                            help = 'Input pickle file for open sea')
    parser.add_argument(   '--coast_fileV3', '-cV3',
                            type = str,
                            required = True,
                            help = 'Input pickle file for coast')

    


    return parser.parse_args()

args = argument()

import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import sys
import numpy as np

fid = open(args.open_sea_fileV3)
LIST = pickle.load(fid)
fid.close()
TIMES,_,_,MODEL_MEAN_v3,SAT___MEAN_v3,_,_ = LIST


fid = open(args.coast_fileV3)
LIST = pickle.load(fid)
fid.close()
model_coast_v3 = LIST[3]

n=len(TIMES)


fid = open(args.open_sea_fileV2)
LIST = pickle.load(fid)
fid.close()
TIMES_v2,_,_,MODEL_MEAN_v2,SAT___MEAN_v2,_,_ = LIST


fid = open(args.coast_fileV2)
LIST = pickle.load(fid)
fid.close()
model_coast_v2 = LIST[3]

TIMES_v2       =TIMES_v2[0:n]
MODEL_MEAN_v2 = MODEL_MEAN_v2[0:n]
SAT___MEAN_v2 = SAT___MEAN_v2[0:n]
model_coast_v2 = model_coast_v2[0:n]



from basins import V2 as OGS
for isub,sub in enumerate(OGS.P):
    print sub.name
    fig, ax = pl.subplots()
    ax.plot(TIMES,SAT___MEAN_v2[:,isub],'og',label=' SAT')
    ax.plot(TIMES,MODEL_MEAN_v2[:,isub],'-k',label=' V2')
    ax.plot(TIMES,MODEL_MEAN_v3[:,isub],'-r',label=' V3')
    ax.plot(TIMES,model_coast_v2[:,isub],':k',label=' V2_coast')
    ax.plot(TIMES,model_coast_v3[:,isub],':r',label=' V3_coast')
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(0.0, 0.6)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30)
    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename=args.outdir+"/"+'chl' + sub.name + ".png"
    fig.set_size_inches((20,6))
    pl.savefig(outfilename)
