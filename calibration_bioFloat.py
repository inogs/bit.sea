# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
#OUTPUTS
#BGC-CLASS4-O2-CORR-PROF-BASIN
#BGC-CLASS4-NIT-CORR-PROF-BASIN
#BGC-CLASS4-CHL-CORR-PROF-BASIN
# of CMEMS-Med-biogeochemistry-ScCP-1.0.pdf

import scipy.io.netcdf as NC
import numpy as np
import os
from commons.time_interval import TimeInterval

from profiler import *
import basins.OGS as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS

M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)


maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

SUBlist=[OGS.wes, OGS.eas]
nSUB = len(SUBlist)
BGC_CLASS4_CHL_CORR_PROF_BASIN = np.zeros((nSUB,), dtype=np.float32)

modelvarname = 'P_i'

for isub, sub in enumerate(SUBlist):
    Profilelist=bio_float.FloatSelector(FLOATVARS[modelvarname],T_INT,sub)
    nP = len(Profilelist)
    
    PROFILEStimeSERIES    = np.zeros((nP,),dtype=[('corr',np.float32),('rmse',np.float32)])
    PROFILEStimeSERIES[:] = np.nan
    
    for ip, p in enumerate(Profilelist):
        singlefloatmatchup = M.getMatchups([p], nav_lev, modelvarname,read_adjusted=True)
        if singlefloatmatchup.number() > 1 :
            PROFILEStimeSERIES['corr'][ip] = singlefloatmatchup.correlation()
            PROFILEStimeSERIES['rmse'][ip] = singlefloatmatchup.RMSE()
    
    BGC_CLASS4_CHL_CORR_PROF_BASIN[isub] = np.nanmean(PROFILEStimeSERIES['corr'])



# These lines just to prove that DOXY and NITRATE do not have ADJUSTED values
import sys
sys.exit()

for modelvarname in ['O2o','N3n']:
    Profilelist=bio_float.FloatSelector(FLOATVARS[modelvarname],T_INT, OGS.med)
    m     = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True )
    m_raw = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=False)
    print modelvarname, "Raw and adjusted values: ", m_raw.number(), m.number()

#Then,
#BGC-CLASS4-O2-CORR-PROF-BASIN
#BGC-CLASS4-NIT-CORR-PROF-BASIN
# are not provided.