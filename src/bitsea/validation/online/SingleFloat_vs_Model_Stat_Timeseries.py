import argparse
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path
from bitsea.utilities.argparse_types import date_from_str
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.
    Produces a file, containing timeseries for some statistics, for each wmo.
    In the outputdir, two new directories will be created, in order to store the output of check.
    - nitrate_check/
    - chla_check/
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
                                required = True)
    
    parser.add_argument(   '--basedir', '-b',
                                type = existing_dir_path,
                                required = True,
                                help = """ PROFILATORE dir, already generated""")    

    parser.add_argument(   '--outdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = "")
    parser.add_argument(   '--date','-d',
                                type = date_from_str,
                                required = True,
                                help = 'start date in yyyymmdd format')

    return parser.parse_args()

args = argument()

import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.instruments import superfloat as bio_float
from bitsea.instruments.matchup_manager import Matchup_Manager
from bitsea.instruments.var_conversions import FLOATVARS
from bitsea.commons.layer import Layer
from bitsea.basins.region import Rectangle
from bitsea.validation.deliverables.metrics2 import find_DCM, find_WBL,find_NITRICL
from bitsea.validation.deliverables.metrics2 import find_OMZ, find_maxO2
from bitsea.validation.deliverables.metrics import find_NITRICL_dz_max
from bitsea.validation.online.metrics import find_WLB
from bitsea.validation.online.SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
from bitsea.basins import V2 as OGS
from datetime import datetime
from dateutil.relativedelta import relativedelta
from bitsea.instruments import check
from bitsea.Float.oxygen_saturation import oxy_sat



BASEDIR = args.basedir
OUTDIR = args.outdir
Check_obj_nitrate = check.check(OUTDIR / "nitrate_check/")
Check_obj_chl     = check.check(OUTDIR / "chla_check/")


TheMask=Mask(args.maskfile, loadtmask=False)

Graphic_DeltaT = relativedelta(months=18)
datestart = args.date -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")

TL = TimeList.fromfilenames(None, BASEDIR / "PROFILES/","ave*.nc")
TI = TimeInterval.fromdatetimes(datestart,args.date)
ALL_PROFILES = bio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))


layer=Layer(0,200)
layer300=Layer(0,350)
layer1000=Layer(200,1000)

VARLIST = ['P_l','N3n','O2o']
Adj = [True,True,False]
extrap = [True,False,True]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','dNit_dz','CM','O2o_sat','OMZ','max_O2']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m
iz300 = TheMask.getDepthIndex(350)+1 # Max Index for depth 300m for Nitracl def
iz10 = TheMask.getDepthIndex(10.8)+1
iz1000 = TheMask.getDepthIndex(1000)+1 # Max Index for depth 1000

for ivar, var_mod in enumerate(VARLIST):
    var = FLOATVARS[var_mod]
    if var_mod == "N3n": Check_obj = Check_obj_nitrate
    if var_mod == "P_l": Check_obj = Check_obj_chl
    if var_mod == "O2o": Check_obj = None
    Profilelist = bio_float.FloatSelector(var, TI, Rectangle(-6,36,30,46))
    wmo_list=bio_float.get_wmo_list(Profilelist)
    for iwmo, wmo in enumerate(wmo_list):
        OUTFILE = OUTDIR / f"{var_mod}_{wmo}.nc"
        print (OUTFILE, flush=True)
        list_float_track=bio_float.filter_by_wmo(Profilelist,wmo)
        nTime = len(list_float_track)
        A_float = np.zeros(( nTime, nStat), np.float32 ) * np.nan
        A_model = np.zeros(( nTime, nStat), np.float32 ) * np.nan

        for itime in range(nTime):
            if "gm200" in locals(): del gm200
            if "gm300" in locals(): del gm300
            p=list_float_track[itime]
            if p.available_params.find(var)<0 : continue

            Pres,Profile,Qc=p.read(var)

            try:

                GM = M.getMatchups2([p], TheMask.zlevels, var_mod, interpolation_on_Float=False,checkobj=Check_obj, extrapolation=extrap[ivar])

            except:
                print (p.ID()  + " not found in " + str(BASEDIR))
                continue


            if GM.number() == 0 :
                print (p.ID() + " excluded")
                continue
            gm200 = GM.subset(layer)
            gm300 = GM.subset(layer300)
            gm1000=GM.subset(layer1000)
 

            nLevels = gm200.number()
            izmax = min(nLevels,iz200)
  
            nLevels300 = gm300.number()
            izmax300 = min(nLevels300,iz300)

            nLevels1000 = gm1000.number()
            izmax1000 = min(nLevels1000,iz1000)

            # INTEGRAL 
            A_float[itime,0] =  np.nansum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
            A_model[itime,0] =  np.nansum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral

            # SURF VALUE
            A_float[itime,5] = gm200.Ref[0] # Surf Value
            A_model[itime,5] = gm200.Model[0] # Surf Value

           # DCM/MWB
            if (var_mod == "P_l"):
                if ( ( p.time.month >= 4. )  & ( p.time.month <= 10 )):
                    A_float[itime,7], A_float[itime,2] = find_DCM(gm200.Ref  ,gm200.Depth) # CM, DCM
                    A_model[itime,7], A_model[itime,2] = find_DCM(gm200.Model,gm200.Depth) # CM, DCM
            
                if (p.time.month in [1,2,3] ):
                    A_float[itime,3] = find_WLB(gm200.Ref  ,gm200.Depth) # WLB
                    A_model[itime,3] = find_WLB(gm200.Model,gm200.Depth) # WLB

           # NITRACL1/NITRACL2 
            if (var_mod == "N3n"):
                # NOTA: level 350
                A_float[itime,4] = find_NITRICL(gm300.Ref  ,gm300.Depth) # Nitricline
                A_model[itime,4] = find_NITRICL(gm300.Model,gm300.Depth) # Nitricline

                A_float[itime,6] = find_NITRICL_dz_max(gm300.Ref  ,gm300.Depth) # dNit/dz
                A_model[itime,6] = find_NITRICL_dz_max(gm300.Model,gm300.Depth) # Nitricline

            if (var_mod == "O2o"):
                A_float[itime,8] = oxy_sat(p)

                if len(gm1000.Ref) > 1:
                    A_float[itime,9] = find_OMZ(gm1000.Ref, gm1000.Depth) # Oxygen Minimum Zone
                    A_model[itime,9] = find_OMZ(gm1000.Model, gm1000.Depth) # Oxygen Minimum Zone 

                    A_float[itime,10] = find_maxO2(gm300.Ref, gm300.Depth) # Oxygen Max depth
                    A_model[itime,10] = find_maxO2(gm300.Model, gm300.Depth) # Oxygen Max depth


            A_float[itime,1] = gm200.correlation() # Correlation
            A_model[itime,1] = gm200.correlation() # Correlation

        dumpfile(OUTFILE,A_float,A_model,METRICS)

    
