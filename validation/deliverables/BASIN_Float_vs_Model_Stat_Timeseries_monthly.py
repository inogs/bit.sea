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
    In the outputdir, two new directories will be created, in order to store the output of check.
    - nitrate_check/
    - chla_check/

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
from instruments import superfloat as bio_float
from instruments.matchup_manager import Matchup_Manager
from instruments.var_conversions import FLOATVARS
from commons.utils import addsep
from commons.layer import Layer
from profiler_floats import ALL_PROFILES,TL,BASEDIR
from metrics2 import *
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
from basins.V2 import NRT3 as OGS
import commons.timerequestors as requestors
from instruments import check
from float_OXY_saturation import *

OUTDIR = addsep(args.outdir)
Check_obj_nitrate = check.check(OUTDIR + "/nitrate_check/")
Check_obj_chl     = check.check(OUTDIR + "chla_check/")
Check_obj_PhytoC  = check.check(OUTDIR + "/Phyto_C/")
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)
layer300=Layer(0,350)
layer1000=Layer(200,1000)

VARLIST = ['P_l','N3n','O2o','P_c']
Adj = [True,True,True,True]
extrap = [True,False,False,False]
nVar = len(VARLIST)
nSub = len(OGS.basin_list)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf']
METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf','dNit_dz','CM','O2o_sat','OMZ','max_O2']
#METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf','OMZ','max_O2']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
MonthlyRequestors=M.TL.getMonthlist()
nTime = len(MonthlyRequestors)


iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m

A_float = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan
A_model = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan

for ivar, var_mod in enumerate(VARLIST):
# if (ivar == 2):
    var = FLOATVARS[var_mod]
    adj=Adj[ivar]

    if var_mod == "N3n": Check_obj = Check_obj_nitrate
    if var_mod == "P_l": Check_obj = Check_obj_chl
    if var_mod == "O2o": Check_obj = None
    if var_mod == "P_c": Check_obj = Check_obj_PhytoC

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
#                Pres,Profile,Qc=p.read(var,read_adjusted=adj)
                if (var_mod=="P_c"):
                    Pres,Profile,Qc=p.read(var,var_mod="P_c")
                else:
                    Pres,Profile,Qc=p.read(var) #,True)

		if len(Pres) < 10 : continue
#		GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, interpolation_on_Float=False)
                GM = M.getMatchups2([p], TheMask.zlevels, var_mod, interpolation_on_Float=False,checkobj=Check_obj, extrapolation=extrap[ivar])


                if GM.number() == 0 :
                    print p.ID() + " excluded"
                    continue

		gm200 = GM.subset(layer)
                gm300 = GM.subset(layer300)
                gm1000=GM.subset(layer1000)

		nLevels = gm200.number()
		izmax = min(nLevels,iz200)

		Flo[ip,0] = np.nansum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
		Flo[ip,1] = gm200.correlation()
		Flo[ip,5] = gm200.Ref[0] # Surf Value

                Mod[ip,0] = np.sum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
                Mod[ip,1] = gm200.correlation()
                Mod[ip,5] = gm200.Model[0] # Surf Value


		if (VARLIST[ivar] == "P_l"):
		    Flo[ip,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
		    Flo[ip,3] = find_WBL(gm200.Ref  ,gm200.Depth) # WBL
		    Mod[ip,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM
                    Mod[ip,3] = find_WBL(gm200.Model,gm200.Depth) # WBL


		if (VARLIST[ivar] == "N3n"):
                    Flo[ip,4] = find_NITRICL(gm200.Ref  ,gm200.Depth) # Nitricline
		    Mod[ip,4] = find_NITRICL(gm200.Model,gm200.Depth)
  
                if (VARLIST[ivar] == "O2o"):
                    if ( len(gm1000.Model) > 2):
                        Flo[ip,7] = find_OMZ(gm1000.Ref, gm1000.Depth) # Oxygen Minimum Zone
                        Mod[ip,7] = find_OMZ(gm1000.Model, gm1000.Depth) # Oxygen Minimum Zone 
                    else:
                        Flo[ip,7] = np.nan
                        Mod[ip,7] = np.nan

                    Flo[ip,8] = find_maxO2(gm300.Ref, gm300.Depth) # Oxygen Max depth
                    Mod[ip,8] = find_maxO2(gm300.Model, gm300.Depth) # Oxygen Max depth


	    for iStat, sStat in enumerate(METRICS):
		if (iStat == 6): continue
	        A_float[ivar,itime,iSub,iStat] = np.nanmean(Flo[ip,iStat])
		A_model[ivar,itime,iSub,iStat] = np.nanmean(Mod[ip,iStat])

    print var
np.save(OUTDIR + 'Basin_Statistics_FLOAT',A_float)
np.save(OUTDIR + 'Basin_Statistics_MODEL',A_model)
