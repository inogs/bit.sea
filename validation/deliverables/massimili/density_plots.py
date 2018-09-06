import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates density plots
    for matchups with static nutrients dataset
    profiler_RA.py defines paths
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['N1p', 'N3n', 'O2o', 'P_l'] )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = '',
                                help = 'Output Images directory')

    parser.add_argument(   '--minvalue', '-m',
                                type = float,
                                required = False,
                                help = 'Minimum value for plot axes')
    parser.add_argument(   '--maskfile', '-M',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')

    return parser.parse_args()


args = argument()


import os
from profiler import Matchup_Manager, MassimiliReader, ALL_PROFILES,TL, BASEDIR,T_INT
from commons.mask import Mask
import basins.V2 as OGS
from instruments.var_conversions import MASSIMILIVARS as NUTRVARS
from commons.utils import addsep


M = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
N=MassimiliReader()

TheMask = Mask(args.maskfile)
nav_lev = TheMask.zlevels

OUTPUTDIR=addsep(args.outdir)
modelvarname = args.varname
os.system('mkdir -p ' + OUTPUTDIR)

UNITS_DICT={'N1p' : 'mmol P/m$^3$', 
         'N3n' : 'mmol N/m$^3$',
         'O2o' :'mmol O$_2$/m$^3$',
         'P_l' :'mmol CHL/m$^3$'
         }


#for sub in [OGS.alb, OGS.nwm, OGS.lev, OGS.ion]:
for sub in OGS.NRT3:
    print sub.name
    Profilelist=N.Selector(NUTRVARS[modelvarname],T_INT,sub)
    nProfiles=len(Profilelist)
    print nProfiles
    if nProfiles==0: continue

    Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
    fig,ax =  Matchup_basin.densityplot2(modelname='RAN',refname='REF',units=UNITS_DICT[modelvarname],sub=sub.name.upper())
    maxval=max(Matchup_basin.Ref.max(),Matchup_basin.Model.max())

    if args.minvalue is None:
        minval=min(Matchup_basin.Ref.min(),Matchup_basin.Model.min())
    else:
        minval=args.minvalue
    
    ax.set_xlim([minval,maxval])
    ax.set_ylim([minval,maxval])
    ax.set_ylabel('RAN data [' + UNITS_DICT[modelvarname] + ']').set_fontsize(14)
    ax.set_xlabel('REF data [' + UNITS_DICT[modelvarname] + ']').set_fontsize(14)
    ax.set_title(sub.name.upper() + ' - TOT n. matchups= ' + str(Matchup_basin.number()))
    ax.grid(True)
    fig.savefig(OUTPUTDIR+'densplot_'+modelvarname+'_'+sub.name+'.png')
