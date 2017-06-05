import os,sys
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
import numpy as np
from commons.mask import Mask
from basins.region import Region, Rectangle
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
import scipy.io.netcdf as NC
import pylab as pl
from commons.layer import Layer
from profiler import *
from metrics import *
from ncwriter import dumpfile
TheMask=Mask('/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc')
layer=Layer(0,200)

VARLIST = ['P_l','N3n','O2o']
Adj = [True,True,False]
nVar = len(VARLIST)

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
OUTDIR = "/pico/scratch/userexternal/lfeudale/validation/work/output/"
wmo_list=bio_float.get_wmo_list(ALL_PROFILES)

izmax = TheMask.getDepthIndex(200) # Max Index for depth 200m

for wmo in wmo_list:
      OUTFILE = OUTDIR + wmo + ".nc"
      print OUTFILE
      list_float_track=bio_float.filter_by_wmo(ALL_PROFILES,wmo)
      nTime = len(list_float_track)
      A_float = np.zeros((nVar, nTime, nStat), np.float32 )
      A_model = np.zeros((nVar, nTime, nStat), np.float32 )


      for ivar, var_mod in enumerate(VARLIST):
	  var = LOVFLOATVARS[var_mod]
          adj=Adj[ivar]
          for itime in range(nTime):

	      p=list_float_track[itime]
              if p.available_params.find(var)<0 : continue
              Pres,Profile,Qc=p.read(var,read_adjusted=adj)
              if len(Pres) < 1 : continue
              GM = M.getMatchups([p], TheMask.zlevels, var_mod, read_adjusted=adj, interpolation_on_Float=False)  
              gm200 = GM.subset(layer)
              #gm.Model
              #gm.Ref
              #gm.bias()
              
              import sys
#              sys.exit()
	      A_float[ivar,itime,0] =  np.sum(gm200.Ref  *TheMask.dz[:izmax+1])/TheMask.dz[:izmax+1].sum() # Integral
	      A_model[ivar,itime,0] =  np.sum(gm200.Model*TheMask.dz[:izmax+1])/TheMask.dz[:izmax+1].sum() # Integral
		
	      A_float[ivar,itime,1] = gm200.correlation() # Correlation
              A_model[ivar,itime,1] = gm200.correlation() # Correlation

	      if (VARLIST[ivar] == "P_l"):	
                  A_float[ivar,itime,2] = find_DCM(gm200.Ref  ,gm200.Depth)[1] # DCM
                  A_model[ivar,itime,2] = find_DCM(gm200.Model,gm200.Depth)[1] # DCM 

                  A_float[ivar,itime,3] = find_MLD(gm200.Ref  ,gm200.Depth) # MLD 
                  A_model[ivar,itime,3] = find_MLD(gm200.Model,gm200.Depth) # MLD

	      if (VARLIST[ivar] == "N3n"):
                  A_float[ivar,itime,4] = find_NITRICL(gm200.Ref  ,gm200.Depth) # Nitricline
                  A_model[ivar,itime,4] = find_NITRICL(gm200.Model,gm200.Depth) # Nitricline


      dumpfile(OUTFILE,A_float,A_model,VARLIST,METRICS)

##################
#          Pres,Prof,Qc=p.read('CHLA',read_adjusted=True)
#          ii = Pres<=200 ;
#
#	  print p.time.strftime("%Y%m%d")
