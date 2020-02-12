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

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )

    parser.add_argument(   '--inputfile', '-i',
                            type = str,
                            required = True,
                            default = 'export_data_ScMYValidation_plan.pkl',
                            help = 'Input pickle file')

    parser.add_argument(   '--runname', '-r',
                            type = str,
                            required = False,
                            default = 'RAN',
                            help = 'Name of the run for legend')

    parser.add_argument(   '--dep', '-d',
                            type = str,
                            required =False,
                            default = "",
                            help = ''' Depth of layer for chl average'''
                            )

    return parser.parse_args()

args = argument()

fid = open(args.inputfile)
LIST = pickle.load(fid)
fid.close()
TIMES                          = LIST[0]
BGC_CLASS4_CHL_RMS_SURF_BASIN  = LIST[1]
BGC_CLASS4_CHL_BIAS_SURF_BASIN = LIST[2]
MODEL_MEAN                     = LIST[3]
SAT___MEAN                     = LIST[4]
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = LIST[5]
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= LIST[6]


from basins import OGS
for isub,sub in enumerate(OGS.P):
    print sub.name
    fig, ax = pl.subplots()
    ax.plot(TIMES,SAT___MEAN[:,isub],'og',label=' SAT')
    ax.plot(TIMES,MODEL_MEAN[:,isub],'-k',label=args.runname)
    ax.set_ylabel(sub.name.upper() + ' DepAve ' + args.dep + 'm' +  ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(0.0, 0.6)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30)
    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename=args.outdir+"/"+'chl' + sub.name + args.dep + ".png"
    pl.savefig(outfilename)
    #sys.exit()
