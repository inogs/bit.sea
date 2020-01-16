
import numpy as np
import matplotlib.pyplot as pl
import argparse
import pickle
import matplotlib.dates as mdates
from basins import V2
from profileruns_2015 import runList,colorList,DICTlegend
#
def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required =True,
                            default = "DA_FLOAT_SAT/Winter/",
                            help = ''' Input dir of .pkl files produced by ScMYvalidation'''
                            )

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')

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

    parser.add_argument(   '--startdate', '-s',
                            type = str,
                            required = True,
                            help = 'start date for plot')

    parser.add_argument(   '--endtdate', '-e',
                            type = str,
                            required = True,
                            help = 'end date for plot')


    return parser.parse_args()

args = argument()

TIMES = {}
MODEL_MEAN = {}
MODEL_MEANallp = {}
SAT___MEAN = {}
for run in runList:
    filenameopen = args.inputdir + '/RUN_' + run + '/' + args.open_sea_file
    fid = open(filenameopen)
    LIST = pickle.load(fid)
    fid.close()
    TIMES[run],_,_, \
    MODEL_MEAN[run],SAT___MEAN[run],_,_, \
    MODEL_MEANallp[run] = LIST



from basins import V2 as OGS
for isub,sub in enumerate(OGS.P):
    print sub.name
    fig, ax = pl.subplots(figsize=[12,5.5])
#
    overall_isub = (OGS.P.basin_list).index(sub)
#    fig,ax = pl.subplots()
    ax.plot(TIMES[run],SAT___MEAN[run][:,isub],'og',label=' SAT')
    for run in runList:
        ax.plot(TIMES[run],MODEL_MEAN[run][:,isub], \
        '-',color=colorList[run],label=DICTlegend[run])
        #ax.plot(TIMES[run],MODEL_MEANallp[run][:,isub], \
        #'-',color=colorList[run],label=run)
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(0.0, 0.6)
    #ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=20)
    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename=args.outdir+"/"+'chlruns' + sub.name + ".png"
    pl.savefig(outfilename)
