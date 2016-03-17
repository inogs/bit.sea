import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates Vertical profiles
    profiler_RA.py defines paths
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['N1p', 'N3n', 'O2o'] )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = '',
                                help = 'Output Images directory')


    return parser.parse_args()


args = argument()

import pylab as pl
import scipy.io.netcdf as NC
import numpy as np
import os
import sys
from commons.time_interval import TimeInterval

from profiler_RA import *

import basins.OGS as OGS
from instruments.instruments import static_Selector
from commons.layer import Layer
M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)

maskfile    = os.getenv("MASKFILE");
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

OUTPUTDIR=args.outdir
modelvarname = args.varname
os.system('mkdir -p ' + OUTPUTDIR)

UNITS_DICT={'N1p' : 'mmol P/m$^3$', 
         'N3n' : 'mmol N/m$^3$',
         'O2o' :'mmol O$_2$/m$^3$' 
         }


for sub in [OGS.alb, OGS.nwm, OGS.lev, OGS.ion]:
#for sub in OGS.P:
    print sub.name
    Profilelist=static_Selector(modelvarname,T_INT,sub)
    Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)


    fig=pl.figure(num=None,dpi=100,facecolor='w',edgecolor='k')
    ax=fig.add_subplot(1,1,1)
    pl.plot(Matchup_basin.Ref,Matchup_basin.Depth,'.r',label='REF')
    pl.plot(Matchup_basin.Model,Matchup_basin.Depth,'.k',label='RAN')
    pl.gca().invert_yaxis()
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    ax.set_xlabel('[' + UNITS_DICT[modelvarname]+ ']').set_fontsize(14)

    ax.set_ylabel('[m]').set_fontsize(14)
    ax.set_title(sub.name.upper() + ' ' + modelvarname)
    fig.savefig(OUTPUTDIR+'vertprof_'+modelvarname+'_'+sub.name+'.png')    

    print sub.name,'correlation= ',Matchup_basin.correlation()




