import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates Vertical profiles
    for matchups with static nutrients dataset
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
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')

    return parser.parse_args()


args = argument()

import matplotlib.pyplot as pl
from commons.mask import Mask
import os
from profiler_RA import *
import basins.OGS as OGS
from instruments.var_conversions import NUTRVARS
from static.Nutrients_reader import NutrientsReader

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
N=NutrientsReader()

TheMask = Mask(args.maskfile)
nav_lev = TheMask.zlevels
OUTPUTDIR=args.outdir
modelvarname = args.varname
os.system('mkdir -p ' + OUTPUTDIR)

UNITS_DICT={'N1p' : 'mmol P/m$^3$', 
         'N3n' : 'mmol N/m$^3$',
         'O2o' :'mmol O$_2$/m$^3$' 
         }


for sub in OGS.P: # do profiles for the sub-basins (needed also for the table)
    Profilelist=N.Selector(NUTRVARS[modelvarname],T_INT,sub)
    nProfiles=len(Profilelist)
    if nProfiles == 0 :
        # no figure generation
        print sub.name,'correlation= ',"nan with 0 profiles and 0 matchups"
        continue
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

    print sub.name,'correlation= %f with %d profiles and %d matchups' %(Matchup_basin.correlation(), nProfiles, Matchup_basin.number() )




