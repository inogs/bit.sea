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


OUTDIR = '/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/OUT_MATCHUP_' + RUN + '/'
MASKFILE = '/pico/scratch/userexternal/ateruzzi/DA_COAST_17/wrkdir/MODEL/meshmask.nc'
TheMask = Mask(MASKFILE)
nav_lev = TheMask.zlevels

area = 'tyr1'

VARLIST = ['N1p','N3n','P_l']
VARLIST = ['N3n']

#SEASLIST = ['winter','spring','summer','autumn']
SEASLIST = ['win-spr','sum-aut']

DICT_seas = {
    'winter': TimeInterval('20130101','20130331','%Y%m%d'),
    'spring': TimeInterval('20130401','20130630','%Y%m%d'),
    'summer': TimeInterval('20130701','20130930','%Y%m%d'),
    'autumn': TimeInterval('20131001','20131231','%Y%m%d'),
    'win-spr': TimeInterval('20130101','20130630','%Y%m%d'),
    'sum-aut': TimeInterval('20130701','20131231','%Y%m%d'),
    }


N = IspraReader()
layer = Layer(0,3)
dictArea = {
    'adriatic' : V2.ComposedBasin('adr',[V2.adr1,V2.adr2]), 
    'adr1'     : V2.adr1, 
    'adr2'     : V2.adr2, 
    'tyr1'     : V2.tyr1, 
    'tyr2'     : V2.tyr2, 
    'ion3'     : V2.ion3
           } 



nMetrics = 6
STATS = np.ones((len(VARLIST),len(SEASLIST),nMetrics))*np.nan
print(RUN)
print(area)
for ivar,var in enumerate(VARLIST):
    print('----------------')
    print(var)
    varname=var_conversions.ISPRAVARS[var]
    for iseas,seas in enumerate(SEASLIST):
        print('................')
        print(seas)
        TI = DICT_seas[seas]
        Profilelist = N.Selector(varname,TI,dictArea[area])
        for p in Profilelist:
        #    if (p.name()=='12-T2_21') or (p.name()=='12-T2_24'): #N1p
        #        print('removing station 12-T2_21 or 24')
        #        Profilelist.remove(p)
            if (p.name()=='12-T2_21') or (p.name()=='12-T2_24'): #N3n
                print('removing station 12-T2_21 or 24')
                Profilelist.remove(p)
        #    if (p.name()=='12-T2_21A') or (p.name()=='12-T2_24A') or (p.name()=='12-T2_24') or (p.name()=='12-T2_21'): #P_l
        #        print('removing station 12-T2_21A or 24A or 24 or 21')
        #        Profilelist.remove(p)
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
        print('len')
        print(len(Llayer.Model))
        print('RMSE')
        print(Llayer.RMSE())
        print('bias')
        print(Llayer.bias())

outfile = OUTDIR + 'stats'+ RUN + '.npy'
np.save(outfile,STATS)

