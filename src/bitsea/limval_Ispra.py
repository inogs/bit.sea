import scipy.io.netcdf as NC
import numpy as np
import os

from commons.layer import Layer
from commons.mask import Mask

from profilerIspra import *
import basins.OGS as OGS
import basins.V2 as V2
from instruments import lovbio_float
from instruments import var_conversions
from static.Ispra_reader import *
from loadIspraDistances import statExcludeInt,statExcludeDist

M = Matchup_Manager(ISPRA_PROFILES,TL,BASEDIR)


OUTDIR = '/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/OUT_MATCHUP_' + RUN + '/'
MASKFILE = '/pico/scratch/userexternal/ateruzzi/DA_COAST_18/wrkdir/MODEL/meshmask.nc'
TheMask = Mask(MASKFILE)
nav_lev = TheMask.zlevels


VARLIST = ['N1p','N3n','P_l']
#VARLIST = ['N3n']

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


AreasList = ['adr1','adr2','tyr1','tyr2','ion3','adriatic','allmed']
dictArea = {
    'adriatic' : V2.ComposedBasin('adr',[V2.adr1,V2.adr2]), 
    'adr1'     : V2.adr1, 
    'adr2'     : V2.adr2, 
    'tyr1'     : V2.tyr1, 
    'tyr2'     : V2.tyr2, 
    'ion3'     : V2.ion3,
    'allmed'   : V2.med
           } 

if __name__ == '__main__':
    N = IspraReader()
    layer = Layer(0,3)

    statExclude = ['12-T2_21A','12-T2_24A','12-T2_24','12-T2_21']

    LIMITS = np.ones((len(VARLIST),len(SEASLIST),len(AreasList)))*np.nan
    print(RUN)
    for iarea,area in enumerate(AreasList):
        print('================')
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
                print(len(Profilelist))
                iexc = 0
                iexcInt = 0
                iexcDist = 0
                newList = []
                for p in Profilelist:
                    #if p.name() in statExclude:
                    #    print('removing station 12-T2_21A or 24A or 24 or 21')
                    #    Profilelist.remove(p)
                    #    iexc += 1
                    condexc1 = p.name().startswith('12-T2_21')
                    condexc4 = p.name().startswith('12-T2_24')
                    condexcI = p.name() in statExcludeInt
                    condexcD = p.name() in statExcludeDist
                    condexc = condexc1 | condexc4 | condexcI | condexcD
                    
                    if condexc1:
                        print(p.name())
                        #Profilelist.remove(p)
                        iexc += 1
                    if condexc:
                        continue
                    else:
                        newList.append(p)
                #L = M.getMatchups(Profilelist, nav_lev, var)
                L = M.getMatchups(newList, nav_lev, var)
                Llayer = L.subset(layer)
                LIMITS[ivar,iseas,iarea] = 1.2*Llayer.Model.max()
                print(len(Profilelist),iexc,iexcInt,iexcDist)
                print(len(newList))
                print(LIMITS[ivar,iseas,iarea])
    
    outfile = OUTDIR + 'limits'+ RUN + '.npy'
    np.save(outfile,LIMITS)

