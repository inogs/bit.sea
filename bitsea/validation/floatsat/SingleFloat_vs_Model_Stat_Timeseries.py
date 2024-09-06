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
from profiler_floatsat import ALL_PROFILES,TL,BASEDIR
from metrics import *
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
import basins.V2 as OGS

OUTDIR = addsep(args.outdir)
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)

VARLIST = ['P_l']
Adj = [True,True,False]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)

iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m

for wmo in wmo_list:
    OUTFILE = OUTDIR + wmo + ".nc"
    print OUTFILE
    list_float_track=bio_float.filter_by_wmo(ALL_PROFILES,wmo)
    nTime = len(list_float_track)
    A_float = np.zeros((nVar, nTime, nStat), np.float32 ) * np.nan
    A_model = np.zeros((nVar, nTime, nStat), np.float32 ) * np.nan


    for ivar, var_mod in enumerate(VARLIST):
        var = LOVFLOATVARS[var_mod]
        adj=Adj[ivar]
        for itime in range(nTime):
            p=list_float_track[itime]
            if p.available_params.find(var)<0 : continue
            Pres,Profile,Qc=p.read(var,read_adjusted=adj)
            if len(Pres) < 10 : continue
            GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, \
                               interpolation_on_Float=False)
            gm200 = GM.subset(layer)
            nLevels = gm200.number()
            izmax = min(nLevels,iz200)
 
            A_float[ivar,itime,0] =  np.sum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
            A_model[ivar,itime,0] =  np.sum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral

            A_float[ivar,itime,1] = gm200.correlation() # Correlation
            A_model[ivar,itime,1] = gm200.correlation() # Correlation

	    A_float[ivar,itime,5] = gm200.Ref[0] # Surf Value
            A_model[ivar,itime,5] = gm200.Model[0] # Surf Value


            if (VARLIST[ivar] == "P_l"):
                A_float[ivar,itime,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
                A_model[ivar,itime,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM
                DCM2_float = find_DCM2(gm200.Ref  ,gm200.Depth)[1] # DCM
                DCM2_model = find_DCM2(gm200.Model,gm200.Depth)[1] # DCM
#               print  "FLOAT MAX: " , A_float[ivar,itime,2]
#		print  "MODELMAX: " , A_model[ivar,itime,2]

                A_float[ivar,itime,3] = find_MLD(gm200.Ref  ,gm200.Depth) # MLD
                A_model[ivar,itime,3] = find_MLD(gm200.Model,gm200.Depth) # MLD

            if (VARLIST[ivar] == "N3n"):
                A_float[ivar,itime,4] = find_NITRICL(gm200.Ref  ,gm200.Depth) # Nitricline
                A_model[ivar,itime,4] = find_NITRICL(gm200.Model,gm200.Depth) # Nitricline
#	    import sys
#	    sys.exit()

    dumpfile(OUTFILE,A_float,A_model,VARLIST,METRICS)
