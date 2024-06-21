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
                                choices = ['P_l','kd490','P1l','P2l','P3l','P4l','RRS412','RRS443','RRS490','RRS510','RRS555','RRS670'],
                                help = ''' model var name'''
                                )
    parser.add_argument(   '--zone', '-z',
                                type = str,
                                required =False,
                                default = "Med",
                                help = ''' Areas to generate the STATISTICS mean. std, bias and RMSD with respect satellite: Med or rivers'''
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

if (args.zone == "Med"):
    from basins import V2 as OGS
if (args.zone == "rivers"):
    from basins import RiverBoxes as OGS

TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD,CORR,NUMB = LIST

model_label=' MODEL'

if (args.var =="kd490"):
    units="[m$^{-1}$]"

elif (args.var.startswith('RRS')):
    units="[st$^{-1}$]"

else:
    units="[mg/m$^3$]"

var_label = SAT_VARS[args.var] + " " + units

vmin=0.0
vmax=0.25
color="tab:green"
lightcolor="palegreen"

if (args.var == "P_l"): 
    vmin=0.0
    if (args.zone == "Med"):
        vmax=0.6
    if (args.zone == "rivers"):
        vmax=1.0 

if (args.var == "kd490"):
    vmin=0.02
    vmax=0.09
    color="tab:blue"
    lightcolor="lightsteelblue"

if (args.var.startswith('RRS')):
    vmin=0.0
    vmax=0.025


for isub,sub in enumerate(OGS.P):
    if (sub.name == 'atl') : continue
    print (sub.name)
    fig, ax = pl.subplots()
    fig.set_size_inches(12,4)
    ax.plot(TIMES,SAT___MEAN[:,isub],'o',label=' SAT',color=color)
    ax.fill_between(TIMES,SAT___MEAN[:,isub]-SAT____STD[:,isub],SAT___MEAN[:,isub]+SAT____STD[:,isub],color=lightcolor)
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
    if (sub.name=="Po"):
       vmax=4.0
    else:
       vmax=1.0
    pl.ylim(vmin,vmax)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=12)
    outfilename="%s%s_%s_STD.png"  %(OUTDIR, args.var, sub.name)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=12)
    pl.tight_layout()
    pl.savefig(outfilename)
