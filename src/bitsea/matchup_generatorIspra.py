import scipy.io.netcdf as NC
import numpy as np
import basins.V2 as V2

from commons.layer import Layer
from commons.mask import Mask

from profilerIspra import *
from instruments import var_conversions
from static.Ispra_reader import *
from limval_Ispra import AreasList as AreasListLim
from limval_Ispra import VARLIST as VARLISTlim
from limval_Ispra import SEASLIST as SEASLISTlim
from loadIspraDistances import statExcludeInt,statExcludeDist


M = Matchup_Manager(ISPRA_PROFILES,TL,BASEDIR)


OUTDIR = '/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/OUT_MATCHUP_' + RUN + '/'
MASKFILE = '/pico/scratch/userexternal/ateruzzi/DA_COAST_18/wrkdir/MODEL/meshmask.nc'
TheMask = Mask(MASKFILE)
nav_lev = TheMask.zlevels
LIMFILE = '/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/OUT_MATCHUP_18/limits18.npy'
LIMITS = np.load(LIMFILE)

area = 'tyr1'

VARLIST = ['N1p','N3n','P_l']
#VARLIST = ['N3n']

#DictLimit = {
#    'N1p': 0.3,
#    'N3n': 10,
#    'P_l': 3
#            }

#SEASLIST = ['winter','spring','summer','autumn']
SEASLIST = ['win-spr','sum-aut']
#SEASLIST = ['win-spr']

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
    'ion3'     : V2.ion3,
    'allmed'   : V2.med
           } 

indarea = AreasListLim.index(area)


# some stations to be excluded (not found..?) --> TO BE RESOLVED
statExclude = ['12-T2_21A','12-T2_24A','12-T2_24','12-T2_21']


nMetrics = 6
STATS = np.ones((len(VARLIST),len(SEASLIST),nMetrics))*np.nan
STATSall = np.ones((len(VARLIST),len(SEASLIST),nMetrics))*np.nan
print(RUN)
print(area)
for ivar,var in enumerate(VARLIST):
    print('----------------')
    print(var)
    varname = var_conversions.ISPRAVARS[var]
    indvar = VARLISTlim.index(var)
    #varlimit = DictLimit[var]
    for iseas,seas in enumerate(SEASLIST):
        indseas = SEASLISTlim.index(seas)
        varlimit = LIMITS[indvar,indseas,indarea]
        print('................')
        print(seas)
        TI = DICT_seas[seas]
        Profilelist = N.Selector(varname,TI,dictArea[area])
        newList = []
        for p in Profilelist:
            condexc1 = p.name().startswith('12-T2_21')
            condexc4 = p.name().startswith('12-T2_24')
            condexcI = p.name() in statExcludeInt
            condexcD = p.name() in statExcludeDist
            condexc = condexc1 | condexc4 | condexcI | condexcD
            if condexc:
                continue
            else:
                newList.append(p)
        #L = M.getMatchups(Profilelist, nav_lev, var)
        L = M.getMatchups(newList, nav_lev, var)
        Llayer = L.subset(layer)
        #allvalue
        outfile = OUTDIR + 'modvals_all' + RUN + seas + var + area + '.npy'
        np.save(outfile,Llayer.Model)
        outfile = OUTDIR + 'refvals_all' + seas + var + area + '.npy'
        np.save(outfile,Llayer.Ref)
        outfile = OUTDIR + 'latitudes_all' + seas + var + area + '.npy'
        np.save(outfile,Llayer.Lat)
        outfile = OUTDIR + 'longitudes_all' + seas + var + area + '.npy'
        np.save(outfile,Llayer.Lon)
        outfile = OUTDIR + 'names_all' + seas + var + area + '.npy'
        np.save(outfile,Llayer.name)
        STATSall[ivar,iseas,0] = len(Llayer.Model)
        STATSall[ivar,iseas,1] = Llayer.Model.std()
        STATSall[ivar,iseas,2] = Llayer.Ref.std()
        STATSall[ivar,iseas,3] = Llayer.correlation()
        STATSall[ivar,iseas,4] = Llayer.RMSE()
        STATSall[ivar,iseas,5] = Llayer.bias()
        print('len all')
        print(len(Llayer.Model))
        print('RMSE all')
        print(Llayer.RMSE())
        print('bias all')
        print(Llayer.bias())
        # limvalue
        Lvalmax = Llayer.limmaxref(varlimit)
        Lvalue  = Lvalmax.limminmodel(0.001)
        outfile = OUTDIR + 'modvals_' + RUN + seas + var + area + '.npy'
        np.save(outfile,Lvalue.Model)
        outfile = OUTDIR + 'refvals_' + seas + var + area + '.npy'
        np.save(outfile,Lvalue.Ref)
        outfile = OUTDIR + 'latitudes_' + seas + var + area + '.npy'
        np.save(outfile,Lvalue.Lat)
        outfile = OUTDIR + 'longitudes_' + seas + var + area + '.npy'
        np.save(outfile,Lvalue.Lon)
        outfile = OUTDIR + 'names_' + seas + var + area + '.npy'
        np.save(outfile,Lvalue.name)
        STATS[ivar,iseas,0] = len(Lvalue.Model)
        STATS[ivar,iseas,1] = Lvalue.Model.std()
        STATS[ivar,iseas,2] = Lvalue.Ref.std()
        STATS[ivar,iseas,3] = Lvalue.correlation()
        STATS[ivar,iseas,4] = Lvalue.RMSE()
        STATS[ivar,iseas,5] = Lvalue.bias()
        print('len')
        print(len(Lvalue.Model))
        print('RMSE')
        print(Lvalue.RMSE())
        print('bias')
        print(Lvalue.bias())

outfile = OUTDIR + 'stats_all'+ RUN + area + '.npy'
np.save(outfile,STATSall)
outfile = OUTDIR + 'stats'+ RUN + area + '.npy'
np.save(outfile,STATS)

