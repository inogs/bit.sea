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

M = Matchup_Manager(ISPRA_PROFILES,TL,BASEDIR)


OUTDIR = '/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/DATI_ISPRA/SKILL_DIAG/OUT_MATCHUP_' + RUN + '/'
MASKFILE = '/pico/scratch/userexternal/ateruzzi/DA_COAST_17/wrkdir/MODEL/meshmask.nc'
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
                for p in Profilelist:
                #    if (p.name()=='12-T2_21') or (p.name()=='12-T2_24'): #N1p
                #        print('removing station 12-T2_21 or 24')
                #        Profilelist.remove(p)
                #    if (p.name()=='12-T2_21') or (p.name()=='12-T2_24'): #N3n
                #        print('removing station 12-T2_21 or 24')
                #        Profilelist.remove(p)
                    if (p.name()=='12-T2_21A') or (p.name()=='12-T2_24A') or (p.name()=='12-T2_24') or (p.name()=='12-T2_21'): #P_l
                        print('removing station 12-T2_21A or 24A or 24 or 21')
                        Profilelist.remove(p)
                L = M.getMatchups(Profilelist, nav_lev, var)
                Llayer = L.subset(layer)
                LIMITS[ivar,iseas,iarea] = 1.2*Llayer.Model.max()
                print(LIMITS[ivar,iseas,iarea])
    
    outfile = OUTDIR + 'limits'+ RUN + '.npy'
    np.save(outfile,LIMITS)

