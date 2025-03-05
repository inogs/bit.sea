import argparse
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path
def argument():
    parser = argparse.ArgumentParser(description = '''

    Generates two files ( model and ref)
    Basin_Statistics_FLOAT.npy
    Basin_Statistics_MODEL.npy
    containing [(nVar, nTime, nSub, nStat)] arrays.
    The metrics are:
    These arrays will be used in the next step to generate the following metrics:

    CHL-PROF-D-CLASS4-PROF-CORR-BASIN
    NIT-PROF-D-CLASS4-PROF-CORR-BASIN
     DO-PROF-D-CLASS4-PROF-CORR-BASIN 

    The following step will be in 
    BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py
    and the results will be displayed in FIG.IV.6 and TABLE IV.3

    In the outputdir, 3 new directories will be created, in order to store the output of check.
    - nitrate_check/
    - chla_check/
    - Phyto_C/

    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
                                required = True,
                                help = ''' Path of maskfile''')
    parser.add_argument(   '--basedir', '-b',
                                type = existing_dir_path,
                                required = True,
                                help = """ PROFILATORE dir, already generated""")
    parser.add_argument(   '--outdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
from bitsea.commons.mask import Mask
from bitsea.instruments import superfloat as bio_float
from bitsea.instruments.matchup_manager import Matchup_Manager
from bitsea.instruments.var_conversions import FLOATVARS
from bitsea.commons.utils import nanmean_without_warnings
from bitsea.commons.layer import Layer
from bitsea.validation.deliverables.metrics import find_DCM, find_WBL,find_NITRICL
from bitsea.validation.deliverables.metrics import find_OMZ, find_maxO2
from bitsea.basins.V2 import NRT3 as OGS
from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.instruments import check
import datetime
from bitsea.basins.region import Rectangle
from bitsea.commons import timerequestors

OUTDIR = args.outdir
BASEDIR = args.basedir

TL = TimeList.fromfilenames(None, BASEDIR / "PROFILES/","ave*.nc")
deltaT= datetime.timedelta(hours=12)
TI = TimeInterval.fromdatetimes(TL.Timelist[0] - deltaT, TL.Timelist[-1] + deltaT)
ALL_PROFILES = bio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


Check_obj_nitrate = check.check(OUTDIR / "nitrate_check")
Check_obj_chl     = check.check(OUTDIR / "chla_check")
Check_obj_PhytoC  = check.check(OUTDIR / "Phyto_C")
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)
layer300=Layer(0,350)
layer1000=Layer(200,1000)

class variable():
    def __init__(self, name, extrap, check_obj):
        ''' Arguments:
        * name *  string, like N3n
        * extrap * logical
        * check_obj * a check object defined in instruments.check
        '''
        self.name = name
        self.extrap = extrap
        self.check_obj = check_obj


P_l = variable('P_l', True, Check_obj_chl)
N3n = variable('N3n', False, Check_obj_nitrate)
O2o = variable('O2o', False, None)
P_c = variable('P_c', False, Check_obj_PhytoC )
VARLIST = [P_l,N3n,O2o,P_c]


nVar = len(VARLIST)
nSub = len(OGS.basin_list)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf']
METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','nProf','dNit_dz','CM','O2o_sat','OMZ','max_O2']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
MonthlyRequestors=M.TL.getMonthlist()
nTime = len(MonthlyRequestors)


iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m

A_float = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan
A_model = np.zeros((nVar, nTime, nSub, nStat), np.float32 ) * np.nan

for ivar, V in enumerate(VARLIST):
    var_mod = V.name
    var = FLOATVARS[var_mod]

    if var_mod == "N3n": Check_obj = Check_obj_nitrate
    if var_mod == "P_l": Check_obj = Check_obj_chl
    if var_mod == "O2o": Check_obj = None
    if var_mod == "P_c": Check_obj = Check_obj_PhytoC

    for itime, Req in enumerate(MonthlyRequestors):
        if Req.time_interval.end_time > TL.timeinterval.end_time :
            Req.time_interval.end_time = TL.timeinterval.end_time
        print (Req)
        for iSub, Sub in enumerate(OGS.basin_list):
            BASIN_PROFILES_float_raw = bio_float.FloatSelector(var,Req.time_interval,Sub)
            BASIN_PROFILES_float = bio_float.remove_bad_sensors(BASIN_PROFILES_float_raw,var)

            A_float[ivar,itime,iSub,6] = len(BASIN_PROFILES_float)
            A_model[ivar,itime,iSub,6] = len(BASIN_PROFILES_float)
            if len(BASIN_PROFILES_float) == 0: continue
    
            Flo = np.zeros((len(BASIN_PROFILES_float), nStat), np.float32 ) * np.nan
            Mod = np.zeros((len(BASIN_PROFILES_float), nStat), np.float32 ) * np.nan
            for ip, p in enumerate(BASIN_PROFILES_float):
                if p.available_params.find(var)<0 : continue
                if (var_mod=="P_c"):
                    Pres,Profile,Qc=p.read(var,var_mod="P_c")
                    print (p.ID())
                else:
                    Pres,Profile,Qc=p.read(var) #,True)

                if len(Pres) < 10 : continue
                GM = M.getMatchups2([p], TheMask.zlevels, var_mod, interpolation_on_Float=False,checkobj=V.check_obj, extrapolation=V.extrap)

                if GM.number() == 0 :
                    print (p.ID() + " excluded")
                    continue

                gm200 = GM.subset(layer)
                gm300 = GM.subset(layer300)
                gm1000=GM.subset(layer1000)

                nLevels = gm200.number()
                izmax = min(nLevels,iz200)

                Flo[ip,0] = np.nansum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
                if gm200.Ref.std()>0 : Flo[ip,1] = gm200.correlation()

                Flo[ip,5] = gm200.Ref[0] # Surf Value


                Mod[ip,0] = np.nansum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
                if gm200.Ref.std()>0 : Mod[ip,1] = gm200.correlation()

                Mod[ip,5] = gm200.Model[0] # Surf Value


                if (var_mod == "P_l"):
                    if ( ( Req.month >= 4. )  & ( Req.month <= 10 )):
                        Flo[ip,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
                        Mod[ip,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM
                    if (Req.month in [1,2,3] ):
                        Flo[ip,3] = find_WBL(gm200.Ref  ,gm200.Depth) # WBL
                        Mod[ip,3] = find_WBL(gm200.Model,gm200.Depth) # WBL


                if (var_mod == "N3n"):
                    Flo[ip,4] = find_NITRICL(gm200.Ref  ,gm200.Depth) # Nitricline
                    Mod[ip,4] = find_NITRICL(gm200.Model,gm200.Depth)
  
                if (var_mod == "O2o"):
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
#                A_float[ivar,itime,iSub,iStat] = np.nanmean(Flo[ip,iStat])
#                A_model[ivar,itime,iSub,iStat] = np.nanmean(Mod[ip,iStat])
                A_float[ivar,itime,iSub,iStat] = nanmean_without_warnings(Flo[:,iStat])
                A_model[ivar,itime,iSub,iStat] = nanmean_without_warnings(Mod[:,iStat])
#                print (A_float[3,:,6,6])

outfile=OUTDIR / "Basin_Statistics_FLOAT.npy"
print(outfile)
np.save(outfile,A_float)
outfile=OUTDIR / "Basin_Statistics_MODEL.npy"
print(outfile)
np.save(outfile,A_model)
