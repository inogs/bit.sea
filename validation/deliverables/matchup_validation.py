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
    parser.add_argument(   '--area', '-a',
                                type = str,
                                default = 'OpenSea',
                                required = False,
                                help = "Choices of the area: OpenSea or Coast")
    return parser.parse_args()

args = argument()

import numpy as np
#from profiler_RA_N import NutrientsReader, Matchup_Manager, ALL_PROFILES, TL, BASEDIR, T_INT
from profiler_RA import Matchup_Manager, ALL_PROFILES, TL, BASEDIR, T_INT
from commons.mask import Mask
import basins.V2 as OGS
from static.climatology import DatasetInfo
from static.Nutrients_reader import NutrientsReader
from static.Carbon_reader import CarbonReader
#from instruments.var_conversions import NUTRVARS as NUTRVARS
#from instruments.var_conversions import CARBONVARS 
from commons.utils import addsep
from commons.layer import Layer

from commons.utils import writetable

M = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
N=NutrientsReader()
C=CarbonReader()

TheMask = Mask(args.maskfile)
nav_lev = TheMask.zlevels
OUTDIR = addsep(args.outdir)
area = args.area

# Define coastal area:
mask200_2D = TheMask.mask_at_level(200.0)
mask0_2D = TheMask.mask_at_level(0.0)
coastmask=mask0_2D & (~mask200_2D)

#LayerList = [Layer(0,1000), Layer(0,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300)]
#LayerList = [Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000), Layer(0,1000)]
LayerList = [Layer(0,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000), Layer(0,1000)]
LayerList_Coast = [Layer(0,60), Layer(60,200), Layer(0,200)]

#SUBLIST=[OGS.nwm, OGS.tyr, OGS.lev, OGS.ion]
SUBLIST = OGS.P.basin_list
SUBLIST.remove(SUBLIST[-1])

if ( area=='Coast' ): LayerList=LayerList_Coast

print LayerList

nSub, nLayer = len(SUBLIST), len(LayerList)
rows_names=[sub.name for sub in SUBLIST]
column_names = [layer.string() for layer in LayerList]


#for modelvarname in ["N1p","N3n","N4n","N5s","O2o","DIC","ALK","pH","pCO2"]:
# Choose here Nutrients or Carbonate System:
for modelvarname in ["N1p","N3n","N4n","N5s","O2o"]:
#for modelvarname in ["DIC","ALK","pH","pCO2"]:
    var, Dataset = DatasetInfo(modelvarname)
    fname = OUTDIR + modelvarname
    STAT = np.zeros((nSub,nLayer,7),np.float32)*np.nan
    NUMB = np.zeros((nSub,nLayer),np.int32)
    NPROF = np.zeros((nSub,1),np.int32)
    
    for isub, sub in enumerate(SUBLIST):
# INITIALIZE THE LIST OF POINTS PER AREA
        Profilelist_Coast = []   # LIST of points in COASTAL AREA
        Profilelist_OpenSea = [] # LIST of points in OPEN SEA area

        Profilelist_all=Dataset.Selector(var,T_INT,sub)
        nProfiles=len(Profilelist_all)
        print sub.name, nProfiles

# select OPEN SEA and COASTAL AREA:
        Lon = np.zeros((nProfiles,), np.float64)*np.nan
        Lat = np.zeros((nProfiles,), np.float64)*np.nan
        for ip, p in enumerate(Profilelist_all):
#            rec_sum=rec_sum + len(p.profile)
            Lon[ip] = p.lon
            Lat[ip] = p.lat

            ix,iy=TheMask.convert_lon_lat_to_indices(Lon[ip],Lat[ip])
            if (coastmask[iy,ix] == True):
               print "Last depth in m is: "    
               print p.pres[-1]
               Profilelist_Coast.append(p) # point in COASTAL AREA
            else:
               Profilelist_OpenSea.append(p) # point in OPEN SEA



        if ( area=='OpenSea' ): Profilelist=Profilelist_OpenSea
        else:
             Profilelist=Profilelist_Coast
             LayerList=LayerList_Coast

        print "AREA: " + area
        nProfiles=len(Profilelist)
        print sub.name, nProfiles

###############

        NPROF[isub,0]=nProfiles
        Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
        for ilayer, layer in enumerate(LayerList):
            m_layer = Matchup_basin.subset(layer)
            print "m_layer.Model"
            print m_layer.Model
            npoints = m_layer.number()
            NUMB[isub,ilayer] = npoints
            if npoints>3:
                print "npoints=", npoints
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



        
