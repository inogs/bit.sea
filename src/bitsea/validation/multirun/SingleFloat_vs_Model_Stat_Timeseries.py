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
from profiler_2015 import ALL_PROFILES,TL,BASEDIR
from metrics import *
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
import basins.V2 as OGS

OUTDIR = addsep(args.outdir)
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)
layer300=Layer(0,300)

VARLIST = ['Chla','N3n','O2o']
VARLIST = ['P_l','N3n','O2o']
VARLIST = ['P_l','N3n']
#VARLIST = ['N3n']
Adj = [True,True,False]
Adj = [True,True]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal', \
    'dNit_dz','Nit_prof','Nit_m10','MLD3','meanNcl']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)

iz200 = TheMask.getDepthIndex(200)+1 # Max Index for depth 200m
iz10 = TheMask.getDepthIndex(10.8)
iz10 = 9
iz10 = TheMask.getDepthIndex(50)

iz10true = TheMask.getDepthIndex(10)

print wmo_list
for wmo in wmo_list:
#   if (wmo=="6901769"):
    print '---------------------------'
    print wmo
    OUTFILE = OUTDIR + wmo + ".nc"
    print OUTFILE
    list_float_track=bio_float.filter_by_wmo(ALL_PROFILES,wmo)
    nTime = len(list_float_track)
    A_float = np.zeros((nVar, nTime, nStat), np.float32 ) * np.nan
    A_model = np.zeros((nVar, nTime, nStat), np.float32 ) * np.nan
#    import sys
#    sys.exit()

    for ivar, var_mod in enumerate(VARLIST):
        print var_mod
    #  if ( ivar == 1 ):
        var = LOVFLOATVARS[var_mod]
        adj=Adj[ivar]
        for itime in range(nTime):
        # for itime in range(1):
            p=list_float_track[itime]
            if p.available_params.find(var)<0 : continue
            Pres,Profile,Qc=p.read(var,read_adjusted=adj)
            if len(Pres) < 10 : continue
            GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, interpolation_on_Float=False)
            # Filter out the N3n spike in surface:
            if (var_mod == "N3n"):
#                for k, kk in enumerate(GM.Ref): 
#                    if (kk<=0.011): 
#                       GM.Ref[k]=np.nan
                # if (GM.Ref[0]>0.44): GM.Ref[:]=np.nan
                if (GM.Ref[0]>2): continue
            if (var_mod == "P_l") or (var_mod == 'Chla'):
                # if (GM.Ref[0]>0.45): GM.Ref[:]=np.nan
                if (GM.Ref[0]>0.45): continue

            gm200 = GM.subset(layer)
            nLevels = gm200.number()
            izmax = min(nLevels,iz200)
            if (var_mod=='P_l') | (var_mod=='Chla'):
                gm300 = GM.subset(layer300)
            if (var_mod == "N3n"):
                layer600=Layer(0,600)
                gm600 = GM.subset(layer600)

            A_float[ivar,itime,0] =  np.nansum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
            A_model[ivar,itime,0] =  np.nansum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral

            A_float[ivar,itime,1] = gm200.correlation() # Correlation
            A_model[ivar,itime,1] = gm200.correlation() # Correlation
            
            A_float[ivar,itime,5] = gm200.Ref[0] # Surf Value
            A_model[ivar,itime,5] = gm200.Model[0] # Surf Value

            A_float[ivar,itime,8] = np.nansum(gm200.Ref[:iz10true]  *TheMask.dz[:iz10true])/TheMask.dz[:iz10true].sum() # Integral
            A_model[ivar,itime,8] = np.nansum(gm200.Model[:iz10true]*TheMask.dz[:iz10true])/TheMask.dz[:iz10true].sum() # Integral

            #if (var_mod == "P_l") or (var_mod == 'Chla'):
            if (var_mod == "N3n"):
                A_float[ivar,itime,5] = gm200.Ref[iz10] # Value at 10m
                A_model[ivar,itime,5] = gm200.Model[iz10] # Value at 10m

            #print "------------"
            #print iz10
            #print gm200.Ref[0]
            #print gm200.Ref[10]
            #print gm200.Ref[iz10]
            #print "------------"


            if ( ( p.time.month >= 4. )  & ( p.time.month <= 10 )):
#                print "DCM " + str(p.time.month)
                A_float[ivar,itime,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
                A_model[ivar,itime,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM
                #DCM2_float = find_DCM2(gm200.Ref  ,gm200.Depth)[1] # DCM
                #DCM2_model = find_DCM2(gm200.Model,gm200.Depth)[1] # DCM
#               print  "FLOAT MAX: " , A_float[ivar,itime,2]
#		print  "MODELMAX: " , A_model[ivar,itime,2]
            if ((( p.time.month >=1 ) & ( p.time.month <= 3 )) | ( p.time.month == 1 )):
#                print "MLB " + str(p.time.month)
                A_float[ivar,itime,3] = find_MLD(gm200.Ref  ,gm200.Depth) # MLD
                A_model[ivar,itime,3] = find_MLD(gm200.Model,gm200.Depth) # MLD
                A_float[ivar,itime,9] = find_MLD(gm300.Ref  ,gm300.Depth) # MLD3
                A_model[ivar,itime,9] = find_MLD(gm300.Model,gm300.Depth) # MLD3

            if (var_mod == "N3n"):
                A_float[ivar,itime,4],_ = find_NITRICL(gm200.Ref  ,gm200.Depth) # Nitricline
                A_model[ivar,itime,4],_ = find_NITRICL(gm200.Model,gm200.Depth) # Nitricline

                A_float[ivar,itime,6] = find_NITRICL_dz_max(gm200.Ref  ,gm200.Depth) # dNit/dz
                A_model[ivar,itime,6] = find_NITRICL_dz_max(gm200.Model,gm200.Depth) # Nitricline

                A_float[ivar,itime,7],izNcl_ref = find_NITRICL(gm600.Ref  ,gm600.Depth) # Nitricline searched also at great depth
                A_model[ivar,itime,7],izNcl_model = find_NITRICL(gm600.Model,gm600.Depth) # Nitricline

                A_float[ivar,itime,10] = np.nansum(gm600.Ref[:izNcl_ref]  *TheMask.dz[:izNcl_ref])/TheMask.dz[:izNcl_ref].sum() # Integral
                A_model[ivar,itime,10] = np.nansum(gm600.Model[:izNcl_model]*TheMask.dz[:izNcl_model])/TheMask.dz[:izNcl_model].sum() # Integral
#                print "NITRICL"


    dumpfile(OUTFILE,A_float,A_model,VARLIST,METRICS)
    import sys
#    sys.exit()
