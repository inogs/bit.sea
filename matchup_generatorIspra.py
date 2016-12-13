import scipy.io.netcdf as NC
import numpy as np
import os

from commons.layer import Layer
from commons.mask import Mask

from profiler import *
import basins.OGS as OGS
import basins.V2 as V2
from instruments import lovbio_float
from instruments import var_conversions
from static.Ispra_reader import *

M = Matchup_Manager(ISPRA_PROFILES,TL,BASEDIR)

RUN = '15'

OUTDIR = '/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/OUT_MATCHUP/'
MASKFILE = '/pico/scratch/userexternal/ateruzzi/DA_COAST_15/wrkdir/MODEL/meshmask.nc'
TheMask = Mask(MASKFILE)
nav_lev = TheMask.zlevels

VARLIST = ['N1p','N3n','P_l']

SEASLIST = ['winter','spring','summer','autumn']

DICT_seas = {
    'winter': TimeInterval('20130101','20130331','%Y%m%d'),
    'spring': TimeInterval('20130401','20130630','%Y%m%d'),
    'summer': TimeInterval('20130701','20130930','%Y%m%d'),
    'autumn': TimeInterval('20131001','20131231','%Y%m%d'),
    }


N = IspraReader()
layer = Layer(0,3)
adriatic = V2.ComposedBasin('adr',[V2.adr1,V2.adr2])


nMetrics = 6
STATS = np.ones((len(VARLIST),len(SEASLIST),nMetrics))*np.nan
for ivar,var in enumerate(VARLIST):
    print('----------------')
    print(var)
    varname=var_conversions.ISPRAVARS[var]
    for iseas,seas in enumerate(SEASLIST):
        print('................')
        print(seas)
        TI = DICT_seas[seas]
        Profilelist = N.Selector(varname,TI,adriatic)
        L = M.getMatchups(Profilelist, nav_lev, var)
        Llayer = L.subset(layer)
        outfile = OUTDIR + 'modvals_' + RUN + seas + var + '.npy'
        np.save(outfile,Llayer.Model)
        outfile = OUTDIR + 'refvals_' + RUN + seas + var + '.npy'
        np.save(outfile,Llayer.Ref)
        STATS[ivar,iseas,0] = len(Llayer.Model)
        STATS[ivar,iseas,1] = Llayer.Model.std()
        STATS[ivar,iseas,2] = Llayer.Ref.std()
        STATS[ivar,iseas,3] = Llayer.correlation()
        STATS[ivar,iseas,4] = Llayer.RMSE()
        STATS[ivar,iseas,5] = Llayer.bias()

outfile = OUTDIR + 'stats'+ RUN + '.npy'
np.save(outfile,STATS)

