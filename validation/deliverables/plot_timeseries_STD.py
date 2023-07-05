import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot timeseries for QUID.
    Reference doc: CMEMS-Med-QUID-006-008-V2-V1.0.docx
    Usually Figure IV.2
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
                            help = 'Input pickle file for open sea')
    parser.add_argument(   '--var', '-v',
                                type = str,
                                required = True,
                                choices = ['P_l','kd490','P1l','P2l','P3l','P4l'],
                                help = ''' model var name'''
                                )
    


    return parser.parse_args()

args = argument()

import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from commons.utils import addsep
from basins import V2 as OGS
from instruments.var_conversions import SAT_VARS

OUTDIR=addsep(args.outdir)
fid = open(args.inputfile,'rb')
LIST = pickle.load(fid)
fid.close()

TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD,CORR = LIST

model_label=' MODEL'

if (args.var =="kd490"):
    units="[m$^{-1}$]"
else:
    units="[mg/m$^3$]"

var_label = SAT_VARS[args.var] + " " + units

vmin=0.0
vmax=0.25

if (args.var == "P_l"): 
    vmin=0.0
    vmax=0.6

if (args.var == "kd490"):
    vmin=0.02
    vmax=0.09

for isub,sub in enumerate(OGS.P):
    if (sub.name == 'atl') : continue
    print (sub.name)
    fig, ax = pl.subplots()
    ax.plot(TIMES,SAT___MEAN[:,isub],'og',label=' SAT')
    ax.fill_between(TIMES,SAT___MEAN[:,isub]-SAT____STD[:,isub],SAT___MEAN[:,isub]+SAT____STD[:,isub],color='palegreen')
    ax.plot(TIMES,MODEL_MEAN[:,isub],'-k',label=model_label)
    ax.plot(TIMES,MODEL_MEAN[:,isub]-MODEL__STD[:,isub],':k')
    ax.plot(TIMES,MODEL_MEAN[:,isub]+MODEL__STD[:,isub],':k')

    ax.set_ylabel("%s - %s" %(sub.name.upper(), var_label  ) ).set_fontsize(14)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    ax.tick_params(axis='both', labelsize=12)
#    pl.ylim(0.0, np.max(MODEL_MEAN[:,isub]+MODEL__STD[:,isub]) * 1.2 )
    pl.ylim(vmin,vmax)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=10)
    outfilename="%s%s_%s_STD.png"  %(OUTDIR, args.var, sub.name)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=10)
    pl.savefig(outfilename)
