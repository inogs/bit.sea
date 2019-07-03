import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    Produces a file, containing timeseries for some statistics, for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
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
from profiler import ALL_PROFILES,TL,BASEDIR
from metrics import *
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
import basins.V2 as OGS

OUTDIR = addsep(args.outdir)
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)
layer300=Layer(0,350)

VARLIST = ['P_l','N3n','O2o']
Adj = [True,True,False]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','dNit_dz']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)

iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m
iz300 = TheMask.getDepthIndex(350)+1 # Max Index for depth 300m for Nitracl def
iz10 = TheMask.getDepthIndex(10.8)+1

class matchup_senza_obs():
    def __init__(self, p, model_varname, max_depth=200):
        '''
        Fake matchup object. It works only with model data.
        '''
        Model_time = M.modeltime(p)
        Modelfile =  M.profilingDir + "PROFILES/" + Model_time.strftime("ave.%Y%m%d-%H:%M:%S.profiles.nc")
        ModelProfile = M.readModelProfile(Modelfile, model_varname, p.ID())
        seaPoints = ~np.isnan(ModelProfile)
        ii =TheMask.zlevels[seaPoints] < max_depth
        self.Depth = TheMask.zlevels[seaPoints][ii]
        self.Model = ModelProfile[seaPoints][ii]
        self.Ref   = np.zeros_like(self.Model)*np.nan

    def number(self):
        return len(self.Model)
    def correlation(self):
        return np.nan


for wmo in wmo_list:
    print wmo_list
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
	    try:
                GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, interpolation_on_Float=False)
             
                gm200 = GM.subset(layer)
                gm300 = GM.subset(layer300)

            except:
                gm200=matchup_senza_obs(p, var_mod)
                gm300=matchup_senza_obs(p, var_mod,max_depth=350)
                print "exception"    

            nLevels = gm200.number()
            izmax = min(nLevels,iz200)
  
            nLevels300 = gm300.number()
            izmax300 = min(nLevels300,iz300)

            # INTEGRAL 
            A_float[ivar,itime,0] =  np.nansum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
            if ( np.isnan(gm200.Ref *TheMask.dz[:izmax]).all() == True ): A_float[ivar,itime,0] =  np.nan
            A_model[ivar,itime,0] =  np.nansum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral

            # SURF VALUE
	    A_float[ivar,itime,5] = gm200.Ref[0] # Surf Value
            if ( VARLIST[ivar] == "N3n" ): print A_float[ivar,itime,5]
            A_model[ivar,itime,5] = gm200.Model[0] # Surf Value

           # DCM/MWB
            if (VARLIST[ivar] == "P_l"):
              if ( ( p.time.month >= 4. )  & ( p.time.month <= 10 )):

                A_float[ivar,itime,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
                A_model[ivar,itime,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM
                DCM2_float = find_DCM2(gm200.Ref  ,gm200.Depth)[1] # DCM
                DCM2_model = find_DCM2(gm200.Model,gm200.Depth)[1] # DCM

              if ((( p.time.month >=1 ) & ( p.time.month <= 3 )) | ( p.time.month == 1 )):

                A_float[ivar,itime,3] = find_WLB(gm200.Ref  ,gm200.Depth) # WLB
                A_model[ivar,itime,3] = find_WLB(gm200.Model,gm200.Depth) # WLB

           # NITRACL1/NITRACL2 
            if (VARLIST[ivar] == "N3n"):
                # NOTA: level 350
                A_float[ivar,itime,4] = find_NITRICL(gm300.Ref  ,gm300.Depth) # Nitricline
                A_model[ivar,itime,4] = find_NITRICL(gm300.Model,gm300.Depth) # Nitricline

                A_float[ivar,itime,6] = find_NITRICL_dz_max(gm300.Ref  ,gm300.Depth) # dNit/dz
                A_model[ivar,itime,6] = find_NITRICL_dz_max(gm300.Model,gm300.Depth) # Nitricline

            A_float[ivar,itime,1] = gm200.correlation() # Correlation
            A_model[ivar,itime,1] = gm200.correlation() # Correlation

    dumpfile(OUTFILE,A_float,A_model,VARLIST,METRICS)
