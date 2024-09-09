import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates density plots
    for matchups with static nutrients dataset
    profiler_RA.py defines paths
    Before running it, the profiler_RA is linked to the correct profiler
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['N1p', 'N3n', 'O2o', 'P_l', 'N4n', 'N5s', 'DIC', 'ALK', 'pH', 'pCO2'] )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = '',
                                help = 'Output Images directory')

#    parser.add_argument(   '--profiler', '-p',
#                                type = str,
#                                required = True,
#                                default = '',
#                                help = 'Profiler for the realtive Dataset')

    parser.add_argument(   '--minvalue', '-m',
                                type = float,
                                required = False,
                                help = 'Minimum value for plot axes')

    parser.add_argument(   '--maskfile', '-M',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')

    parser.add_argument(   '--coastness', '-c',
                                type = str,
                                required = True,
                                help = 'COASTNESS list: everywhere, open_sea, coast, supercoastal')
    return parser.parse_args()


args = argument()


import os
import numpy as np
#import __init__
#from profiler_RA_N import *
#from profiler_N4n import *
#from profiler_RA_C import *
from bitsea.static.climatology import DatasetInfo
from bitsea.commons.mask import Mask
import bitsea.basins.V2 as OGS
import bitsea.basins.COASTAL12nm as C12
from bitsea.commons.utils import addsep
from bitsea.static.Nutrients_reader import NutrientsReader
from bitsea.static.Carbon_reader import CarbonReader

# C12.NAd_coastal_basins
#from bitsea.instruments.var_conversions import NUTRVARS, CARBONVARS , SOCAT_VARS
#from bitsea.static.Socat_reader import SocatReader

N=NutrientsReader()
C=CarbonReader()
#S=SocatReader()

TheMask = Mask(args.maskfile)
nav_lev = TheMask.zlevels

# Define coastal area:
mask200_2D = TheMask.mask_at_level(200.0)
mask0_2D = TheMask.mask_at_level(0.0)
coastmask=mask0_2D & (~mask200_2D)

OUTPUTDIR=addsep(args.outdir)
modelvarname = args.varname
coastness = args.coastness 
#input_profiler = args.profiler
#print input_profiler
os.system('mkdir -p ' + OUTPUTDIR)
#from input_profiler import *
from profiler_RA import *
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

UNITS_DICT={'N1p' : 'mmol P/m$^3$', 
         'N3n' : 'mmol N/m$^3$',
         'O2o' :'mmol O$_2$/m$^3$',
         'P_l' :'mmol CHL/m$^3$',
         'N5s' :'mmol SiO$_2$/m$^3$',
         'N4n' :'mmol NH$_4$/m$^3$',
         'ALK' : 'ALK',
         'DIC' : 'DIC',
          'pH' : 'pH',
        'pCO2' : 'pCO2'
         }


var, Dataset = DatasetInfo(modelvarname)

if (coastness == "supercoastal"):
   SUBBASINS = C12.NAd_coastal_basins
else:
   SUBBASINS = OGS.NRT3

#for sub in [OGS.alb, OGS.nwm, OGS.lev, OGS.ion]:
#for sub in OGS.NRT3:
#for sub in C12.NAd_coastal_basins:
for sub in SUBBASINS:
#  if (sub.name =="nwm"):
    print (sub.name)
    Profilelist_OpenSea = [] # LIST of points in OPEN SEA area
    Profilelist_Coast = [] # LIST of points along the coast
    Profilelist_all=Dataset.Selector(var,T_INT,sub)
    nProfiles=len(Profilelist_all)
    print (nProfiles)
    if nProfiles==0: continue

############
# SELECT ONLY OPEN SEA:
    Lon = np.zeros((nProfiles,), np.float64)*np.nan
    Lat = np.zeros((nProfiles,), np.float64)*np.nan

    # CREATE THE LIST OF TRADITIONAL COAST or OPEN_SEA datasets:
    if (coastness=="coast" or coastness=="open_sea"):
        for ip, p in enumerate(Profilelist_all):
#            rec_sum=rec_sum + len(p.profile)
            Lon[ip] = p.lon
            Lat[ip] = p.lat

            ix,iy=TheMask.convert_lon_lat_to_indices(Lon[ip],Lat[ip])
            if (coastmask[iy,ix] == False):
               Profilelist_OpenSea.append(p) # point in OPEN SEA
            if (coastmask[iy,ix] == True):
               Profilelist_Coast.append(p) # point along the coast

#    Profilelist=Profilelist_OpenSea
    if (coastness == "coast"):    Profilelist=Profilelist_Coast
    if (coastness == "open_sea"): Profilelist=Profilelist_OpenSea
#    Profilelist=Profilelist_Coast
    if (coastness == "supercoastal"): Profilelist=Profilelist_all
    if (coastness == "everywhere"):   Profilelist=Profilelist_all
############
    print ("MODELVARNAME: " + modelvarname) 
    Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
    print ("Number of matchups: " + str(Matchup_basin.number()))
    if (Matchup_basin.number() == 0): continue

    fig,ax =  Matchup_basin.densityplot2(modelname='MOD',refname='REF',units=UNITS_DICT[modelvarname],sub=sub.name.upper())
    maxval=max(Matchup_basin.Ref.max(),Matchup_basin.Model.max())

    if args.minvalue is None:
        minval=min(Matchup_basin.Ref.min(),Matchup_basin.Model.min())
    else:
        minval=args.minvalue
    
    ax.set_xlim([minval,maxval])
    ax.set_ylim([minval,maxval])
    ax.set_ylabel('MOD data [' + UNITS_DICT[modelvarname] + ']').set_fontsize(14)
    ax.set_xlabel('REF data [' + UNITS_DICT[modelvarname] + ']').set_fontsize(14)
    ax.set_title(sub.name.upper() + ' - TOT n. matchups= ' + str(Matchup_basin.number()))
    ax.grid(True)
    fig.savefig(OUTPUTDIR+'densplot_'+modelvarname+'_'+sub.name+'.png')
