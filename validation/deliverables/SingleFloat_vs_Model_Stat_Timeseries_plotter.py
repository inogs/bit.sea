import os,sys
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
import numpy as np
from commons.mask import Mask
from basins.region import Region, Rectangle
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
import pylab as pl
import matplotlib.pyplot as plt
from profiler import *
from ncreader import *

def fig_setup(wmo,Lon,Lat,var):
    from layer_integral import coastline

    fig = plt.figure()
    ax0 = plt.subplot2grid((4, 3), (0, 0), colspan=2)
    ax1 = plt.subplot2grid((4, 3), (0, 2))
    ax2 = plt.subplot2grid((4, 3), (1, 0), colspan=3)
    ax3 = plt.subplot2grid((4, 3), (2, 0), colspan=3)
    ax4 = plt.subplot2grid((4, 3), (3, 0), colspan=3)
    axs = [ax0, ax1, ax2, ax3, ax4]

    fig.set_size_inches(10,15)
    c_lon,c_lat=coastline.get()

#    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo_list[j])
    ax0.plot(c_lon,c_lat,'k')
    ax0.plot(Lon,Lat,'r.')
    ax0.set_title("TRAJECTORY of FLOAT " + wmo , color = 'r')
#    ind_max_sup=plotmat[0,:].argmax()
    
#    print Lon[ind_max_sup],Lat[ind_max_sup]
#    ax0.plot(Lon[ind_max_sup],Lat[ind_max_sup],'g.')
#    ax0.plot(Lon[0],Lat[0],'bx')
    ax0.set_xlim([-10,36])
    ax0.set_ylabel("LAT",color = 'k')
    ax0.set_xlabel("LON",color = 'k')

    return fig, axs


TheMask=Mask('/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc')
INDIR = "/pico/scratch/userexternal/lfeudale/validation/work/output/"
OUTDIR = "/pico/scratch/userexternal/lfeudale/validation/work/output/PNG/"

VARLIST = ['P_l','N3n','O2o']
nVar = len(VARLIST)
METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
wmo_list=bio_float.get_wmo_list(ALL_PROFILES)

izmax = TheMask.getDepthIndex(200) # Max Index for depth 200m

for wmo in wmo_list:
      INPUT_FILE = INDIR + wmo + ".nc"
      print INPUT_FILE
      A = ncreader(INPUT_FILE)
      wmo_track_list = bio_float.filter_by_wmo(ALL_PROFILES,wmo)
      nP = len(wmo_track_list)
      Lon = np.zeros((nP,), np.float64)
      Lat = np.zeros((nP,), np.float64)
      for ip, p in enumerate(wmo_track_list):
          Lon[ip] = p.lon
          Lat[ip] = p.lat
      times = [p.time for p in wmo_track_list]
      for var in VARLIST:
	  OUTFILE = OUTDIR + var + "_" + wmo + ".png"
	  print OUTFILE
	  fig, axes = fig_setup(wmo,Lon,Lat,var)
	  if (var == "P_l"): 
	      model, float =A.plotdata(var,'Int_0-200')
              axes[2].plot(times,float,'r')
              axes[2].plot(times,model,'b')
          fig.savefig(OUTFILE)
