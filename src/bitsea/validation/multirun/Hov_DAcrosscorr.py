import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
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
import os,sys
import scipy.io.netcdf as NC
from bitsea.commons.mask import Mask
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.utils import addsep
from bitsea.instruments import matchup_manager
from bitsea.instruments import lovbio_float as bio_float
from bitsea.instruments.var_conversions import LOVFLOATVARS
from bitsea.layer_integral import coastline
from bitsea.instruments.matchup_manager import Matchup_Manager
from bitsea.basins import OGS
import matplotlib.dates as mdates
from bitsea.mhelpers.pgmean import PLGaussianMean
import numpy.ma as ma
import matplotlib.pyplot as pl
from profiler_comparison2015 import *
from profiler_onlymodel import dep1,dep2
from bitsea.timeseries.plot import *


#meanObj11 = PLGaussianMean(5,1.0)

TheMask = Mask.from_file(args.maskfile)
OUTDIR = addsep(args.outdir)

font_s =  15 
label_s = 15

def readModelProfile(filename,var, wmo):
    ncIN = NC.netcdf_file(filename,'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile,:]
    return Profile

def get_leveldep(TheMask,dep):

    i0 = TheMask.getDepthIndex(dep)
    i1 = i0 + 1
      
    diff0 = dep - TheMask.zlevels[i0]
    diff1 = dep - TheMask.zlevels[i1]
    data = [(i0,diff0),(i1,diff1)]
    ix,_ = min(data, key=lambda t: t[1])

    return ix

T_start = DATESTART
T_end   = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start,'%Y%m%d'))
T_end2num   = mpldates.date2num(datetime.strptime(T_end,'%Y%m%d'))
reg1 = [OGS.med]
reg_sn = ['med']

jpk,_,_ = TheMask.shape
max_depth1 = get_leveldep(TheMask,dep1)
max_depth2 = get_leveldep(TheMask,dep2)

MM = Matchup_Manager(ALL_PROFILES,TL,BASEDIR[RUN_DA])


Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)


print wmo_list

for j,wmo in enumerate(wmo_list):
# for j in range(0,1):
    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    nP = len(list_float_track)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)
      
    plotmat = np.zeros([jpk, len(list_float_track)])
    plotmat[:,:] = np.nan
    timelabel_list = list()

    for ip, p in enumerate(list_float_track):
        Lon[ip] = p.lon
        Lat[ip] = p.lat
        timelabel_list.append(p.time)

        plotmat[:,ip] = 


	print var_mod + " " + np.str(len(timelabel_list))

	fig = pl.figure()
	fig.set_size_inches(15,10)

	ax1 = pl.subplot2grid((3, 3), (0, 0), colspan=2)
	ax2 = pl.subplot2grid((3, 3), (0, 2))
	ax3 = pl.subplot2grid((3, 3), (1, 0), colspan=3)
	ax4 = pl.subplot2grid((3, 3), (2, 0), colspan=3)
	
	xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), \
                depths_part)
        plotmat_innov_m=ma.masked_invalid(plotmat_innov[0:nlev_part,:])
        if (var_mod == 'P_l'):
                quadmesh = ax3.pcolormesh(xs, ys, plotmat_innov_m,shading='flat', \
                        vmin=-0.1,vmax=0.1,cmap="Spectral_r")# default is 'flat'
        if (var_mod == 'N1p'):
                quadmesh = ax3.pcolormesh(xs, ys, plotmat_innov_m,shading='flat', \
                        vmin=-0.05,vmax=0.05,cmap="Spectral_r") #,cmap="jet")# default is 'flat'
        if (var_mod == 'N3n'):
                quadmesh = ax3.pcolormesh(xs, ys, plotmat_innov_m,shading='flat', \
                        vmin=-2,vmax=2,cmap="Spectral_r") #,cmap="jet")# default is 'flat'
	ax3.set_xlim([T_start2num,T_end2num])
	ax3.set_xticklabels([])
	ax3.invert_yaxis()
#	ax3.set_title("BFM model: " + plotvarname[ivar], color = 'b')
	ax3.set_ylabel("depth $[m]$",color = 'k',fontsize=font_s)


	xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), \
                depths[nlev_part:-1])
        plotmat_innov_m=ma.masked_invalid(plotmat_innov[nlev_part:-1,:])
        if (var_mod == 'P_l'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_innov_m,shading='flat', \
                        vmin=-0.1,vmax=0.1,cmap="Spectral_r")# default is 'flat'
        if (var_mod == 'N1p'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_innov_m,shading='flat', \
                        vmin=-0.05,vmax=0.05,cmap="Spectral_r") #,cmap="jet")# default is 'flat'
        if (var_mod == 'N3n'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_innov_m,shading='flat', \
                        vmin=-2,vmax=2,cmap="Spectral_r") #,cmap="jet")# default is 'flat'
	ax4.set_xlim([T_start2num,T_end2num])
	ax4.xaxis_date()
	ax4.invert_yaxis()
#	ax3.set_title("BFM model: " + plotvarname[ivar], color = 'b')
	ax4.set_ylabel("depth $[m]$",color = 'k',fontsize=font_s)



        ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax4.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m")) #("%b%y"))
        for tick in ax3.get_xticklabels():
            tick.set_rotation(45)
        
        ax1.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax2.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax3.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax4.tick_params(axis = 'both', which = 'major', labelsize = label_s)

	colorbar_bottom= ax4.get_position().ymin
	colorbar_extent= ax3.get_position().ymax - colorbar_bottom

  #      cbaxes = fig.add_axes([0.95, 0.01, 0.02, 0.95]) 
	cbaxes = fig.add_axes([0.93, colorbar_bottom , 0.02, colorbar_extent])
	cbar=fig.colorbar(quadmesh, cax = cbaxes, )  
	ticklabs = cbar.ax.get_yticklabels()
	cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)


	c_lon,c_lat=coastline.get()

	ax1.plot(c_lon,c_lat,'k')
	ax1.plot(Lon,Lat,'r.')
	ax1.set_title("TRAJECTORY of FLOAT " + p.name() + " - " + var_mod, color = 'r', fontsize=font_s)
	# ind_max_sup=plotmat[0,:].argmax()
        ax1.plot(Lon[0],Lat[0],'bo')
	ax1.plot(Lon[0],Lat[0],'bx')
	ax1.set_xlim([-5.5,36])
#	ax1.set_ylabel("LAT",color = 'k',fontsize=font_s)
#	ax1.set_xlabel("LON",color = 'k')
	extent=4 #degrees
	ipp = len(list_float_track)
	ax2.plot(c_lon,c_lat,'k')
	ax2.plot(Lon,Lat,'ro')
#	ax2.plot(Lon[ind_max_sup],Lat[ind_max_sup],'go')
	ax2.plot(Lon[0],Lat[0],'bo')
	ax2.plot(Lon[0],Lat[0],'bx',markersize=font_s)
	ax2.set_xlim([np.min(Lon[:ipp]) -extent/2, np.max(Lon[:ipp]) +extent/2])
	ax2.set_ylim([np.min(Lat[:ipp]) -extent/2, np.max(Lat[:ipp]) +extent/2])

	fig.savefig(''.join([OUTDIR,'Hov_innovation_',p.name(),'_',var_mod,'.png']))
	pl.close(fig)
