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
                                choices = ['N1p', 'N3n', 'O2o', 'N5s', 'P_l', 'N4n','ALK','DIC','pH','pCO2'] )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = '',
                                help = 'Output Images directory')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')

    parser.add_argument(   '--basedir', '-b',
                                type = str,
                                required = True,
                                help = """ PROFILATORE dir, already generated""")

    parser.add_argument(   '--coastness', '-c',
                                type = str,
                                required = True,
                                help = 'COASTNESS list: everywhere, open_sea, coast')
    parser.add_argument(   '--zone', '-z',
                                type = str,
                                required =False,
                                default = "Med",
                                help = ''' Areas to generate the STATISTICS mean. std, bias and RMSD with respect satellite: Med or rivers''')

    return parser.parse_args()


args = argument()

import matplotlib.pyplot as pl
from bitsea.commons.mask import Mask
import os
import numpy as np
from bitsea.static.Nutrients_reader import NutrientsReader
from bitsea.static.Carbon_reader import CarbonReader
from bitsea.static.climatology import DatasetInfo
from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.basins.region import Rectangle
import datetime
from bitsea.instruments.matchup_manager import Matchup_Manager

N=NutrientsReader()
C=CarbonReader()

BASEDIR = args.basedir
TL = TimeList.fromfilenames(None, BASEDIR +'/'+ "PROFILES/","ave*.nc")
deltaT= datetime.timedelta(hours=12)
T_INT = TimeInterval.fromdatetimes(TL.Timelist[0] - deltaT, TL.Timelist[-1] + deltaT)
ALL_PROFILES = N.Selector(None,T_INT, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
TheMask = Mask.from_file(args.maskfile)
nav_lev = TheMask.zlevels

OUTPUTDIR=args.outdir
modelvarname = args.varname
os.system('mkdir -p ' + OUTPUTDIR)

# Define coastal area:
mask200_2D = TheMask.mask_at_level(200.0)
mask0_2D = TheMask.mask_at_level(0.0)
coastmask=mask0_2D & (~mask200_2D)
coastness = args.coastness
area     = args.zone

UNITS_DICT={'N1p' : 'mmol P/m$^3$', 
         'N3n' : 'mmol N/m$^3$',
         'P_l' :'mmol CHL/m$^3$',
         'O2o' :'mmol O$_2$/m$^3$' ,
         'N4n' : 'mmol NH4/m$^3$',
         'N5s' : 'mmol SiO2/m$^3$',
         'ALK' : 'ALK',
         'DIC' : 'DIC',
         'pH'  : 'pH',
         'pCO2':'pCO2'
         }

var, Dataset = DatasetInfo(modelvarname)

if (area=="Med"):
    from bitsea.basins import V2 as OGS
    SUBBASINS=OGS.P
if (area=="coast12"):
    import bitsea.basins.COASTAL12nm as OGS
    SUBBASINS= OGS.P
if (area=="rivers"):
    from bitsea.basins import RiverBoxes as OGS
    SUBBASINS=OGS.P


#for sub in OGS.P: # do profiles for the sub-basins (needed also for the table)
for sub in SUBBASINS:
    Profilelist_all=Dataset.Selector(var,T_INT,sub)
    nProfiles=len(Profilelist_all)
    if nProfiles == 0 :
        # no figure generation
        print (sub.name,'correlation= ',"nan with 0 profiles and 0 matchups")
        continue

########################
    Profilelist_OpenSea = [] # LIST of points in OPEN SEA area
    Profilelist_Coast   = [] # LIST of points on the COAST

# select OPEN SEA area:
    Lon = np.zeros((nProfiles,), np.float64)*np.nan
    Lat = np.zeros((nProfiles,), np.float64)*np.nan
    for ip, p in enumerate(Profilelist_all):
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
    if (coastness == "everywhere"):   Profilelist=Profilelist_all

    Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)

    fig=pl.figure(num=None,dpi=100,facecolor='w',edgecolor='k')
    ax=fig.add_subplot(1,1,1)
    pl.plot(Matchup_basin.Ref,Matchup_basin.Depth,'.r',label='REF')
    pl.plot(Matchup_basin.Model,Matchup_basin.Depth,'.k',label='MOD')
    pl.gca().invert_yaxis()
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    ax.set_xlabel('[' + UNITS_DICT[modelvarname]+ ']').set_fontsize(14)

    ax.set_ylabel('[m]').set_fontsize(14)
    ax.set_title(sub.name.upper() + ' ' + modelvarname)
    fig.savefig(OUTPUTDIR+'/vertprof_'+modelvarname+'_'+sub.name+'.png')    

    print (sub.name,'correlation= %f with %d profiles and %d matchups' %(Matchup_basin.correlation(), nProfiles, Matchup_basin.number() ))




