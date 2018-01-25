import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    Produces a file, containing timeseries for some statistics, for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
from commons.mask import Mask
from instruments import lovbio_float as bio_float
from instruments.matchup_manager import Matchup_Manager
from instruments.var_conversions import LOVFLOATVARS
from commons.utils import addsep
from commons.layer import Layer
from profiler_floatsat_DAdates import ALL_PROFILES,TL_DA,BASEDIR
from metrics import *
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
import basins.V2 as OGS

OUTDIR = addsep(args.outdir)
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)
layer150=Layer(0,150)

VARLIST = ['P_l']
Adj = [True,True,False]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','Int_0-150']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL_DA,BASEDIR)

MED_PROFILES = bio_float.FloatSelector(None,TL_DA.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)

iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m
iz150 = TheMask.getDepthIndex(150)+1 # Max Index for depth 150m

nTime = len(M.Coupled_List)

A_float = {}
A_model = {}
for wmo in wmo_list:
    #list_float_track=bio_float.filter_by_wmo(ALL_PROFILES,wmo)
    #nTime = len(list_float_track)
    A_float[wmo] = np.zeros((nVar, nTime, nStat), np.float32 ) * np.nan
    A_model[wmo] = np.zeros((nVar, nTime, nStat), np.float32 ) * np.nan

for itime,couple in enumerate(M.Coupled_List):
    intref = {}
    intmod = {}
    in250r = {}
    in250m = {}
    corr = {}
    ref = {}
    mod = {}
    DCMref = {}
    DCMmod = {}
    DCM2r = {}
    DCM2m = {}
    MLDr = {}
    MLDm = {}
    for wmo in wmo_list:
        intref[wmo] = []
        intmod[wmo] = []
        in250r[wmo] = []
        in250m[wmo] = []
        corr[wmo] = []
        ref[wmo] = []
        mod[wmo] = []
        DCMref[wmo] = []
        DCMmod[wmo] = []
        DCM2r[wmo] = []
        DCM2m[wmo] = []
        MLDr[wmo] = []
        MLDm[wmo] = []
#        NTRr[wmo] = []
#        NTRm[wmo] = []
    indexesFloats = couple[1]
    for ifl in indexesFloats:
        p = M.PROFILE_LIST[ifl]
        if p.name() in wmo_list:
            for ivar, var_mod in enumerate(VARLIST):
                var = LOVFLOATVARS[var_mod]
                adj=Adj[ivar]
                floatname = p.name()
                if p.available_params.find(var)<0 : continue
                Pres,Profile,Qc=p.read(var,read_adjusted=adj)
                if len(Pres) < 10 : continue
                GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, \
                        interpolation_on_Float=False)
                gm200 = GM.subset(layer)
                gm150 = GM.subset(layer150)
                nLevels = gm200.number()
                nLevels150 = gm150.number()
                izmax = min(nLevels,iz200)
                izmax150 = min(nLevels,iz150)

                intref[floatname].append(np.sum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum())# Integral
                intmod[floatname].append(np.sum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum())# Integral

                in250r[floatname].append(np.sum(gm150.Ref  *TheMask.dz[:izmax150])/TheMask.dz[:izmax150].sum())# Integral
                in250m[floatname].append(np.sum(gm150.Model*TheMask.dz[:izmax150])/TheMask.dz[:izmax150].sum())# Integral

                corr[floatname].append(gm200.correlation())# Correlation

                ref[floatname].append(gm200.Ref[0])# Surf Value
                mod[floatname].append(gm200.Model[0])# Surf Value

                if (VARLIST[ivar] == "P_l"):
                    DCMref[floatname].append(find_DCM(gm200.Ref  ,gm200.Depth)[1]) # DCM
                    DCMmod[floatname].append(find_DCM(gm200.Model,gm200.Depth)[1])# DCM
                    DCM2r[floatname].append(find_DCM2(gm200.Ref  ,gm200.Depth)[1]) # DCM
                    DCM2m[floatname].append(find_DCM2(gm200.Model,gm200.Depth)[1]) # DCM

                    MLDr[floatname].append(find_MLD(gm200.Ref  ,gm200.Depth))# MLD
                    MLDm[floatname].append(find_MLD(gm200.Model,gm200.Depth))# MLD

#	    import sys
#	    sys.exit()
    for wmo in wmo_list:
        for ivar, var_mod in enumerate(VARLIST):
            A_float[wmo][ivar,itime,0] = np.nanmean(intref[wmo])
            A_model[wmo][ivar,itime,0] = np.nanmean(intmod[wmo]) 

            A_float[wmo][ivar,itime,6] = np.nanmean(in250r[wmo])
            A_model[wmo][ivar,itime,6] = np.nanmean(in250m[wmo])

            A_float[wmo][ivar,itime,1] = np.nanmean(corr[wmo])   # Correlation
            A_model[wmo][ivar,itime,1] = np.nanmean(corr[wmo])   # Correlation

            A_float[wmo][ivar,itime,5] = np.nanmean(ref[wmo]) # Surf Value
            A_model[wmo][ivar,itime,5] = np.nanmean(mod[wmo]) # Surf Value


            if (VARLIST[ivar] == "P_l"):
                A_float[wmo][ivar,itime,2] = np.nanmean(DCMref[wmo]) # DCM
                A_model[wmo][ivar,itime,2] = np.nanmean(DCMmod[wmo]) # DCM

                A_float[wmo][ivar,itime,3] = np.nanmean(MLDr[wmo]) # MLD
                A_model[wmo][ivar,itime,3] = np.nanmean(MLDm[wmo]) # MLD

#	    import sys
#	    sys.exit()

for wmo in wmo_list:
    OUTFILE = OUTDIR + wmo + ".nc"
    print OUTFILE
    dumpfile(OUTFILE,A_float[wmo],A_model[wmo],VARLIST,METRICS)
