# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates:
    BGC_CLASS4_NIT_RMS_LAYER_BASIN
    BGC_CLASS4_PHOS_RMS_LAYER_BASIN
    BGC_CLASS4_O2_RMS_LAYER_BASIN
    BGC_CLASS4_CHL_RMS_LAYER_BASIN ---> not available
    
    of CMEMS-Med- biogeochemistry-ScMYVP-1.0.pdf
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outfile', '-o',
                            type = str,
                            required = True,
                            default = 'export_data_MYValidation_plan_static.pkl',
                            help = 'Output pickle file')


    return parser.parse_args()

args = argument()

import scipy.io.netcdf as NC
import numpy as np
import os
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

SUBlist = [sub for sub in OGS.P]
SUBlist.extend([OGS.wes, OGS.eas])


LAYERLIST = [Layer(0,10), Layer(10,50), Layer(50,100), Layer(100,150), Layer(150,300),Layer(300,600), Layer(600,1000),
             Layer(0,50)] # questo aggiunto in caso di mancanza di dato

nSUB = len(SUBlist)
nLay = len(LAYERLIST)

BGC_CLASS4_NIT_RMS_LAYER_BASIN  = np.ones((nSUB,nLay), dtype=np.float32)*np.NaN
BGC_CLASS4_PHOS_RMS_LAYER_BASIN = np.ones((nSUB,nLay), dtype=np.float32)*np.NaN
BGC_CLASS4_CHL_RMS_LAYER_BASIN  = np.ones((nSUB,nLay), dtype=np.float32)*np.NaN
BGC_CLASS4_O2_RMS_LAYER_BASIN   = np.ones((nSUB,nLay), dtype=np.float32)*np.NaN

VARLIST = ['N1p', 'N3n','O2o']#'chl']

for modelvarname in VARLIST:

    OUTPUT = np.ones((nSUB,nLay), dtype=np.float32)*np.NaN
    for isub, sub in enumerate(SUBlist):
        Profilelist=static_Selector(modelvarname,T_INT,sub)
        nP = len(Profilelist)
        Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
        for ilayer, layer in enumerate(LAYERLIST):
            Mlayer = Matchup_basin.subset(layer)
            OUTPUT[isub,ilayer] = Mlayer.RMSE()
    if modelvarname == 'N1p': BGC_CLASS4_PHOS_RMS_LAYER_BASIN = OUTPUT.copy()
    if modelvarname == 'N3n': BGC_CLASS4_NIT_RMS_LAYER_BASIN  = OUTPUT.copy()
    if modelvarname == 'O2o': BGC_CLASS4_O2_RMS_LAYER_BASIN   = OUTPUT.copy()

import pickle
LIST=[BGC_CLASS4_PHOS_RMS_LAYER_BASIN,
BGC_CLASS4_NIT_RMS_LAYER_BASIN,
BGC_CLASS4_O2_RMS_LAYER_BASIN,
SUBlist,
LAYERLIST]

fid = open(args.outfile,'wb')
pickle.dump(LIST, fid)
fid.close()

