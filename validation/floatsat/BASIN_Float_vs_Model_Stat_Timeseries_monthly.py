import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    Generates in output directory two files ( model and ref) 
    containing [(nVar, nTime, nSub, nStat)] arrays.
    The metrics are:
    These arrays will be used in the next step to generate the following metrics:

    CHL-PROF-D-CLASS4-PROF-CORR-BASIN
    NIT-PROF-D-CLASS4-PROF-CORR-BASIN
     DO-PROF-D-CLASS4-PROF-CORR-BASIN 

    The following step will be in 
    BASIN_Float_vs_Model_Stat_Timeseries_monthly.py
    and the results will be displayed in FIG.IV.6 and TABLE IV.3

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
from basins.V2 import NRT3 as OGS
import commons.timerequestors as requestors

OUTDIR = addsep(args.outdir)
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)

VARLIST = ['P_l'] #,'N3n','O2o']
Adj = [True,True,False]
nVar = len(VARLIST)
nSub = len(OGS.basin_list)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
MonthlyRequestors=M.TL.getMonthlist()
nTime = len(MonthlyRequestors)


iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m

A_float = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan
A_model = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan

for ivar, var_mod in enumerate(VARLIST):
    var = LOVFLOATVARS[var_mod]
    adj=Adj[ivar]
    for itime, Req in enumerate(MonthlyRequestors):
	if Req.time_interval.end_time > TL.timeinterval.end_time : 
	    Req.time_interval.end_time = TL.timeinterval.end_time
	print Req
	for iSub, Sub in enumerate(OGS.basin_list):
	    BASIN_PROFILES_float_raw = bio_float.FloatSelector(var,Req.time_interval,Sub)
            BASIN_PROFILES_float = bio_float.remove_bad_sensors(BASIN_PROFILES_float_raw,var)
            print "RAW  : " + np.str(len(BASIN_PROFILES_float_raw))
            print "FILT.: " + np.str(len(BASIN_PROFILES_float))
	    A_float[ivar,itime,iSub,6] = len(BASIN_PROFILES_float)
            A_model[ivar,itime,iSub,6] = len(BASIN_PROFILES_float)
	    if len(BASIN_PROFILES_float) == 0: continue
	    
	    Flo = np.zeros((len(BASIN_PROFILES_float), nStat), np.float32 ) * np.nan
	    Mod = np.zeros((len(BASIN_PROFILES_float), nStat), np.float32 ) * np.nan
	    for ip, p in enumerate(BASIN_PROFILES_float):
		if p.available_params.find(var)<0 : continue
                Pres,Profile,Qc=p.read(var,read_adjusted=adj)

		if len(Pres) < 10 : continue
		GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, interpolation_on_Float=False)
		gm200 = GM.subset(layer)
		nLevels = gm200.number()
		izmax = min(nLevels,iz200)

		Flo[ip,0] = np.sum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
		Flo[ip,1] = gm200.correlation()
		Flo[ip,5] = gm200.Ref[0] # Surf Value

                Mod[ip,0] = np.sum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
                Mod[ip,1] = gm200.correlation()
                Mod[ip,5] = gm200.Model[0] # Surf Value


		if (VARLIST[ivar] == "P_l"):
		    Flo[ip,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
		    Flo[ip,3] = find_MLD(gm200.Ref  ,gm200.Depth) # MLB
		    Mod[ip,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM
                    Mod[ip,3] = find_MLD(gm200.Model,gm200.Depth) # MLB


		if (VARLIST[ivar] == "N3n"):
                    Flo[ip,4] = find_NITRICL(gm200.Ref  ,gm200.Depth) # Nitricline
		    Mod[ip,4] = find_NITRICL(gm200.Model,gm200.Depth)


	    for iStat, sStat in enumerate(METRICS):
		if (iStat == 6): continue
	        A_float[ivar,itime,iSub,iStat] = np.nanmean(Flo[ip,iStat])
		A_model[ivar,itime,iSub,iStat] = np.nanmean(Mod[ip,iStat])

    print var
np.save(OUTDIR + 'Basin_Statistics_FLOAT',A_float)
np.save(OUTDIR + 'Basin_Statistics_MODEL',A_model)
