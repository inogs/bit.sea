# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>

#OUTPUTS
#BGC-CLASS4-O2-CORR-PROF-BASIN
#BGC-CLASS4-NIT-CORR-PROF-BASIN
#BGC-CLASS4-CHL-CORR-PROF-BASIN
# of CMEMS-Med-biogeochemistry-ScCP-1.0.pdf

import scipy.io.netcdf as NC
import numpy as np
import os
from commons.time_interval import TimeInterval

from profiler import *
import basins.OGS as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
from commons.layer import Layer
M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)

OUTPUT_DIR = "/pico/home/userexternal/gcossari/COPERNICUS/Carbonatic17/CFR_BIOARGO/"

maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

#SUBlist=[OGS.wes, OGS.eas]
SUBlist=[OGS.alb, OGS.nwm, OGS.sww, OGS.swe, OGS.tyr, OGS.ads, OGS.ion, OGS.lev, OGS.aeg ]
nSUB = len(SUBlist)
layer = Layer(0,200)
BGC_CLASS4_CHL_CORR_PROF_BASIN = np.zeros((nSUB,15), dtype=np.float32)
BGC_CLASS4_CHL_RMS_PROF_BASIN = np.zeros((nSUB,15), dtype=np.float32)
BGC_CLASS4_CHL_NUM_PROF_BASIN = np.zeros((nSUB,15), dtype=np.float32)

modelvarname = 'P_i'

MonthlyRequestors=M.TL.getMonthlist() # create time requestor: monthly

StringMonth=['A','M','J','J14','A','S','O','N','D','J15','F','M','A','M','J']
xlabelSI=   [ 0 , 0 , 0 , 0   ,0  , 0 , 0 , 0 , 0 , 0   , 0 , 0 , 1 , 1 ,1]
ylabelSI=   [ 0 , 0 , 0 , 1   ,0  , 0 , 1 , 0 , 0 , 1   , 0 , 0 , 1 , 0 ,0]

import matplotlib.pyplot as pl

for isub, sub in enumerate(SUBlist):
  fig1=pl.figure(num=None,dpi=300,facecolor='w',edgecolor='k')
  for im, Req in enumerate(MonthlyRequestors):
   if im>2: # salto A, M, J,  
    ax1=fig1.add_subplot(4,3,im-2)

    Profilelist=bio_float.FloatSelector(FLOATVARS[modelvarname],Req.time_interval,sub)
    nP = len(Profilelist)
    if nP>1:    
     PROFILEStimeSERIES    = np.zeros((nP,),dtype=[('corr',np.float32),('rmse',np.float32)])
     PROFILEStimeSERIES[:] = np.nan
    
     for ip, p in enumerate(Profilelist):
        singlefloatmatchup = M.getMatchups([p], nav_lev, modelvarname,read_adjusted=True)
        sfm_top = singlefloatmatchup.subset(layer)
        if sfm_top.number() > 1 :
            PROFILEStimeSERIES['corr'][ip] = sfm_top.correlation()
            PROFILEStimeSERIES['rmse'][ip] = sfm_top.RMSE()

            ax1.plot(sfm_top.Model,-sfm_top.Depth,'-k',linewidth=.5)
            ax1.hold(True)
            ax1.plot(sfm_top.Ref,-sfm_top.Depth,'-b',linewidth=.5)
     pl.ylim(-200, 0)
     pl.xlim(0, 1)
     if im==4:
       pl.title(sub.name.upper()).set_fontsize(12)
     if xlabelSI[im]==1:
       ax1.set_xlabel('CHL [mg/m^3]').set_fontsize(10)
       ax1.ticklabel_format(fontsize=8)
#       TickLabels=ax1.get_xticklabels()
#       ax1.set_xticklabels(TickLabels,fontsize=8)
     else:
       ax1.set_xticklabels('',fontsize=8)
     if ylabelSI[im]==1:
       ax1.set_ylabel('depth [m]').set_fontsize(10)
       ax1.ticklabel_format(fontsize=8)
#       TickLabels=ax1.get_yticklabels()
#       ax1.set_yticklabels(TickLabels,fontsize=8)
     else:
       ax1.set_yticklabels('',fontsize=8)
     ax1.grid(True)
#    ax1.text(0,0,StringMonth[im],horizontalalignment='left',verticalalignment='top',fontsize=12, color='black')
     ax1.text(0.98,-190,StringMonth[im],horizontalalignment='right',verticalalignment='bottom',fontsize=14, color='black')
#        transform=ax.transAxes)

     BGC_CLASS4_CHL_CORR_PROF_BASIN[isub,im] = np.nanmean(PROFILEStimeSERIES['corr'])
     BGC_CLASS4_CHL_RMS_PROF_BASIN[isub,im] = np.nanmean(PROFILEStimeSERIES['rmse'])
     BGC_CLASS4_CHL_NUM_PROF_BASIN[isub,im] = (~np.isnan(PROFILEStimeSERIES['corr'])).sum()



  nomefile=OUTPUT_DIR + 'plot_profile_chl_mod_bioargo_' +sub.name + '.png'
  fig1.savefig(nomefile)
  pl.close()


    
nomefile=OUTPUT_DIR + 'BIOARGO_corr_month_subbas.txt'
np.savetxt(nomefile,BGC_CLASS4_CHL_CORR_PROF_BASIN)

nomefile=OUTPUT_DIR + 'BIOARGO_rms_month_subbas.txt'
np.savetxt(nomefile,BGC_CLASS4_CHL_RMS_PROF_BASIN)

nomefile=OUTPUT_DIR + 'BIOARGO_num_month_subbas.txt'
np.savetxt(nomefile,BGC_CLASS4_CHL_NUM_PROF_BASIN)
# These lines just to prove that DOXY and NITRATE do not have ADJUSTED values
import sys
#sys.exit()

for modelvarname in ['P_i','O2o','N3n']:
    Profilelist=bio_float.FloatSelector(FLOATVARS[modelvarname],T_INT, OGS.med)
    WMO = set()
    for p in Profilelist: WMO.add(p._my_float.wmo)
    m     = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True )
    m_raw = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=False)
    print modelvarname
    print "Raw      values: " , m_raw.number()
    print "adjusted values: " , m.number()
    print "Number of files having a void ADJUSTED variable:" , len(Profilelist)
    print "The associated wmo list is :"
    for wmo in WMO:
        print wmo

#Then,
#BGC-CLASS4-O2-CORR-PROF-BASIN
#BGC-CLASS4-NIT-CORR-PROF-BASIN
# are not provided.
