import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--season', '-s',
                                type = str,
                                default = 'annual',
                                required = False,
                                help = "Choices of the seasons are: annual - ann, winter - win, summer - sum")

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

import numpy as np
from bitsea.commons.mask import Mask
from bitsea.static.climatology import DatasetInfo
from bitsea.static.Nutrients_reader import NutrientsReader
from bitsea.static.Carbon_reader import CarbonReader
from bitsea.commons.utils import addsep
from bitsea.commons.layer import Layer
from bitsea.commons.utils import writetable
from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.basins.region import Rectangle
import datetime
from bitsea.instruments.matchup_manager import Matchup_Manager
import pandas as pd
import warnings
warnings.filterwarnings("ignore")


TheMask = Mask.from_file(args.maskfile)
nav_lev = TheMask.zlevels
OUTDIR = addsep(args.outdir)
area = args.zone
coastness = args.coastness
Seas = args.season

N=NutrientsReader()
C=CarbonReader()

# carol added 2025-02-24
BASEDIR = args.basedir
TL = TimeList.fromfilenames(None, BASEDIR +'/'+ "PROFILES/","ave*.nc")
deltaT= datetime.timedelta(hours=12)
T_INT = TimeInterval.fromdatetimes(TL.Timelist[0] - deltaT, TL.Timelist[-1] + deltaT)
ALL_PROFILES = N.Selector(None,T_INT, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
# carol added 2025-02-24

# Define coastal area:
mask200_2D = TheMask.mask_at_level(200.0)
mask0_2D = TheMask.mask_at_level(0.0)
coastmask=mask0_2D & (~mask200_2D)


if (area=="Med"):
    from bitsea.basins import V2 as OGS
    SUBLIST=OGS.P
    SUBLIST.remove(SUBLIST[-1])
if (area=="coast12"):
    import bitsea.basins.COASTAL12nm as OGS
    SUBLIST= OGS.P.basin_list
if (area=="rivers"):
    from bitsea.basins import RiverBoxes as OGS
    SUBLIST =OGS.P.basin_list

LayerList_Coast =[Layer(0,5)]
LayerList= LayerList_Coast
nSub, nLayer = len(SUBLIST), len(LayerList)
rows_names=[sub.name for sub in SUBLIST]
column_names = [layer.string() for layer in LayerList]

#for modelvarname in ["N1p","N3n","O2o"]:
#for modelvarname in ["N1p"]:
for modelvarname in ["P_l", "N1p","N3n","O2o"]: # ""N1p","N3n","N4n","N5s","O2o"]:
#for modelvarname in ["N3n","O2o"]:
    var, Dataset = DatasetInfo(modelvarname)
    fname = OUTDIR + modelvarname
    STAT = np.zeros((nSub,nLayer,7),np.float32)*np.nan
    NUMB = np.zeros((nSub,nLayer),np.int32)
    NPROF = np.zeros((nSub,1),np.int32)
    
    for isub, sub in enumerate(SUBLIST):
        Profilelist_Coast = []   # LIST of points in COASTAL AREA
        Profilelist_OpenSea = [] # LIST of points in OPEN SEA area
        Profilelist_all=Dataset.Selector(var,T_INT,sub)
        nProfiles=len(Profilelist_all)
        Lon = np.zeros((nProfiles,), np.float64)*np.nan
        Lat = np.zeros((nProfiles,), np.float64)*np.nan
        for ip, p in enumerate(Profilelist_all):
            Lon[ip] = p.lon
            Lat[ip] = p.lat
            ix,iy=TheMask.convert_lon_lat_to_indices(Lon[ip],Lat[ip])
            if (coastmask[iy,ix] == True):
               Profilelist_Coast.append(p) # point in COASTAL AREA
            else:
               Profilelist_OpenSea.append(p) # point in OPEN SEA
        if ( area=='OpenSea' ): Profilelist=Profilelist_OpenSea
        else:
             Profilelist=Profilelist_Coast
             LayerList=LayerList_Coast

        print ("AREA: " + area)
        nProfiles=len(Profilelist)

###############

        NPROF[isub,0]=nProfiles
        Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
        for ilayer, layer in enumerate(LayerList):
            m_layer = Matchup_basin.subset_exclude_layerlimits(layer)
            npoints = m_layer.number()
            NUMB[isub,ilayer] = npoints
            if npoints>3:
                print ("npoints=", npoints)
                STAT[isub,ilayer,0]=m_layer.bias()
                STAT[isub,ilayer,1]=m_layer.RMSE()
                STAT[isub,ilayer,2]=m_layer.correlation()
                STAT[isub,ilayer,3]=np.nanmean(m_layer.Model)
                STAT[isub,ilayer,4]=np.nanmean(m_layer.Ref)
                STAT[isub,ilayer,5]=np.nanstd(m_layer.Model)
                STAT[isub,ilayer,6]=np.nanstd(m_layer.Ref)
    writetable(fname + ".bias.txt", STAT[:,:,0], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".rmse.txt", STAT[:,:,1], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".corr.txt", STAT[:,:,2], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".numb.txt", NUMB       , rows_names, column_names, fmt='%1.0f\t')
    writetable(fname + ".nprofiles.txt", NPROF, rows_names, ["TOT PROFILES"], fmt='%2.0f\t')

    writetable(fname + ".ModMEAN.txt", STAT[:,:,3], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".RefMEAN.txt", STAT[:,:,4], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".ModSTD.txt", STAT[:,:,5], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".RefSTD.txt", STAT[:,:,6], rows_names, column_names, fmt='%5.3f\t')

"""

# SEASON CHOICE:
if (Seas == "win" or Seas == "sum"):
    from bitsea.commons.season import season
    S=season()
    S.setseasons(["0101", "0501", "0601", "1001"], ["winter","spring","summer","fall"])
    from bitsea.commons import timerequestors
    from bitsea.commons.Timelist import TimeInterval, TimeList
    if (Seas == "win" ):
        DATESTART = '20190101'
        DATE__END = '20190501'

        T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
    if (Seas == "sum"):
        DATESTART = '20190601'
        DATE__END = '20191101'

    T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
#    TL=TimeList(TIMES)
    from bitsea.commons.utils import writetable

    iSeas=0 # JAN-APR
    CLIM_REQ=timerequestors.Clim_season(iSeas,S)
    iSeas=2 # JUN-SEP
    CLIM_REQ=timerequestors.Clim_season(iSeas,S)
"""

        
