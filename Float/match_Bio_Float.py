import scipy.io.netcdf as NC
import numpy as np
import os
from postproc.Timelist import *
from matchup.matchup import *
from basins.region import *
import instruments.bio_float as Float_Manager
from shared_data import *

import MatchGenerator 
import basins.OGS as OGS

from commons.time_interval import TimeInterval

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')
M = MatchGenerator.Float_Matchup_Manager(T_INT,INPUTDIR,BASEDIR)

maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

MODELVARS={'DOXY':'O2o', \
           'NITRATE':'N3n',  \
           'CHLA':'P_i'}
           

# R1 = Rectangle(-6,10,30,46)
# R2 = Rectangle(10,36,30,46)
# SUBLIST = [R1,R2]

VARLIST=[ 'DOXY','CHLA']


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


for ivar, ref_varname  in enumerate(VARLIST):
    model_varname = MODELVARS[ref_varname]
    print model_varname
    for isub, sub in enumerate(SUBLIST) :
        print sub
        FLOAT_LIST=Float_Manager.FloatSelector(ref_varname, TI, sub)
        Matchup = M.getMatchups(FLOAT_LIST, nav_lev,model_varname,ref_varname)
    
        for ilayer, layer in enumerate(LAYERLIST) :
            print layer
            Mlayer = Matchup.subset(layer)
            
            n = Mlayer.number()
            if n>0:
                STATS[ivar, isub, ilayer, 0] = n 
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



# Provare lo stesso per le fisiche
