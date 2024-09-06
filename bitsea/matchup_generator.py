import scipy.io.netcdf as NC
import numpy as np
import os

from commons.layer import Layer

from profiler import *
import basins.OGS as OGS
from instruments import lovbio_float
from instruments import var_conversions
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

var='O2o'
varname=var_conversions.LOVFLOATVARS[var]
Profilelist=lovbio_float.FloatSelector(varname,T_INT,OGS.tyr)
L = M.getMatchups(Profilelist, nav_lev, var)



import sys
sys.exit()

VARLIST=[ 'O2o','P_i']


SUBLIST = [sub for sub in OGS.P]
LAYERLIST = [Layer(0,10), Layer(10,50), Layer(50,100), Layer(100,150), Layer(150,300),Layer(300,600), Layer(600,1000)]
TI = TimeInterval("20150903", "20150917", "%Y%m%d")
nVars = len(VARLIST)
nSub  = len(SUBLIST)
nLay  = len(LAYERLIST)

METRICS_NAMES=['number of data values',\
               'mean of product',\
               'mean of reference',\
               'mean squared error',\
               'variance of product',\
               'variance of reference']

nMetrics = len(METRICS_NAMES)

STATS = np.ones((nVars,nSub,nLay,nMetrics),np.float32)*np.nan



for ivar, model_varname  in enumerate(VARLIST):
    print model_varname
    for isub, sub in enumerate(SUBLIST) :
        print sub
        Profilelist=instruments.Selector(model_varname, TI, sub)
        Matchup = M.getMatchups(Profilelist, nav_lev,model_varname)
    
        for ilayer, layer in enumerate(LAYERLIST) :
            print layer
            Mlayer = Matchup.subset(layer)
            
            n = Mlayer.number()
            STATS[ivar, isub, ilayer, 0] = n 
            if n>0:
                STATS[ivar, isub, ilayer, 1] = Mlayer.Model.mean()
                STATS[ivar, isub, ilayer, 2] = Mlayer.Ref.mean()
                STATS[ivar, isub, ilayer, 3] = Mlayer.MSE()
                STATS[ivar, isub, ilayer, 4] = np.median(Mlayer.Model)
                STATS[ivar, isub, ilayer, 5] = np.median(Mlayer.Ref)
         

outfilename = "outfile.nc"
ncOUT = NC.netcdf_file(outfilename,'w')
ncOUT.createDimension('nVars', nVars)
ncOUT.createDimension('nSub', nSub)
ncOUT.createDimension('nLay', nLay)
ncOUT.createDimension("metrics",nMetrics     )

ncvar=ncOUT.createVariable('METRICS', 'f', ('nVars','nSub','nLay','metrics'))
ncvar[:] = STATS

ncOUT.close()
