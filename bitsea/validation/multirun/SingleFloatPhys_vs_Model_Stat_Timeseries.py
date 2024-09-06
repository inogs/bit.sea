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
from profiler_phys import ALL_PROFILES,TL,BASEDIR
from metrics import *
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import dumpfile
import basins.V2 as OGS

OUTDIR = addsep(args.outdir)
TheMask=Mask(args.maskfile, loadtmask=False)
layer=Layer(0,200)
layer300=Layer(0,300)

VARLIST = ['votemper','vosaline']
#VARLIST = ['N3n']
Adj = [False,False]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','Val_m10','Termocl1','Termocl2']
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
    OUTFILE = OUTDIR + wmo + "phys.nc"
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
            # Filter on surface value
            # if (var_mod == "P_l") or (var_mod == 'Chla'):
            #     # if (GM.Ref[0]>0.45): GM.Ref[:]=np.nan
            #     if (GM.Ref[0]>0.45): continue

            gm200 = GM.subset(layer)
            nLevels = gm200.number()
            izmax = min(nLevels,iz200)
            layer600=Layer(0,600)
            gm600 = GM.subset(layer600)

            A_float[ivar,itime,0] =  np.nansum(gm200.Ref  *TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral
            A_model[ivar,itime,0] =  np.nansum(gm200.Model*TheMask.dz[:izmax])/TheMask.dz[:izmax].sum() # Integral

            A_float[ivar,itime,1] = gm200.correlation() # Correlation
            A_model[ivar,itime,1] = gm200.correlation() # Correlation
            
            A_float[ivar,itime,2] = gm200.Ref[iz10true] # Value at 10m
            A_model[ivar,itime,2] = gm200.Model[iz10true] # Value at 10m
            

            if (var_mod == "votemper"):
                A_float[ivar,itime,3],_ = find_THERMOCL(gm600.Ref  ,gm600.Depth,14) # Thermocline 14
                A_model[ivar,itime,3],_ = find_THERMOCL(gm600.Model,gm600.Depth,14) # Thermocline

                A_float[ivar,itime,4],_ = find_THERMOCL(gm600.Ref  ,gm600.Depth,16) # Thermocline 16
                A_model[ivar,itime,4],_ = find_THERMOCL(gm600.Model,gm600.Depth,16) # Thermocline




    dumpfile(OUTFILE,A_float,A_model,VARLIST,METRICS)
    import sys
#    sys.exit()
