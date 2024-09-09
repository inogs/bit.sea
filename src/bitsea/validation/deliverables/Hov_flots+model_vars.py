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

import os,sys
import scipy.io.netcdf as NC
from commons.mask import Mask
from instruments import matchup_manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Region, Rectangle
from layer_integral import coastline
import instruments
from instruments import superfloat as bio_float
from instruments.var_conversions import FLOATVARS

import scipy.io.netcdf as NC
import numpy as np
from commons.utils import addsep
import matplotlib.pyplot as pl
from basins import OGS
from validation.online.profileplotter import figure_generator, ncwriter #, add_metadata
import matplotlib
from mhelpers.pgmean import PLGaussianMean
from profiler import *
from instruments.matchup_manager import Matchup_Manager
import matplotlib.dates as mdates

from commons.mask import Mask
from timeseries.plot import *
meanObj11 = PLGaussianMean(5,1.0)
import matplotlib.pyplot as plt
import numpy.ma as ma

TheMask=Mask(args.maskfile)
OUTDIR = addsep(args.outdir)

font_s =  15 
label_s = 15

def my_Hovmoeller_diagram(plotmat, xs,ys, fig=None, ax=None):
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
    quadmesh = ax.pcolormesh(xs, ys, plotmat,shading='gouraud')# default is 'flat'
    #Inform matplotlib that the x axis is made by dates
    ax.xaxis_date()
    ax.invert_yaxis()
    return fig, ax, quadmesh

def readModelProfile(filename,var, wmo):
    ncIN = NC.netcdf_file(filename,'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile,:]
    return Profile

def get_level300(TheMask):

    i0 = TheMask.getDepthIndex(300)
    i1 = i0 + 1
      
    diff0 = 300 - TheMask.zlevels[i0]
    diff1 = 300 - TheMask.zlevels[i1]
    data = [(i0,diff0),(i1,diff1)]
    ix,datamin = min(data, key=lambda t: t[1])

    return ix

T_start = DATESTART
T_end   = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start,'%Y%m%d'))
T_end2num   = mpldates.date2num(datetime.strptime(T_end,'%Y%m%d'))
reg1 = [OGS.med]
reg_sn = ['med']

#max_depth = 26 
#max_depth = 44 # 303m in the 142 levels model
max_depth = get_level300(TheMask)

MM = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
#varname = ['CHLA','DOXY','NITRATE','TEMP','PSAL']
plotvarname = [r'Chl $[mg/m^3]$',r'Oxy $[mmol/m^3]$',r'Nitr $[mmol/m^3]$'] #,r'Temp $[^\circ C]$','Sal']
#read_adjusted = [True,False,False,False,False]
#mapgraph = [3,4,5,1,2]
VARLIST = ['P_l','O2o','N3n']
Adj = [True,False,True]
nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11,1.0)

Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

print wmo_list

for j in range(0,len(wmo_list)):
#    if (j==0): 
      list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo_list[j])
      nP = len(list_float_track)
      Lon = np.zeros((nP,), np.float64)
      Lat = np.zeros((nP,), np.float64)
      NewPres_5m=np.linspace(0,300,61)
      depths=NewPres_5m

      for ivar, var_mod in enumerate(VARLIST):
#       if (ivar==0):
        plotmat = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_model = np.zeros([len(depths), len(list_float_track)])*np.nan
	timelabel_list = list()
        var = FLOATVARS[var_mod]
        adj=Adj[ivar]
        
	for ip, p in enumerate(list_float_track):
#         if ( ip == 91 ):
	    Lon[ip] = p.lon
	    Lat[ip] = p.lat
	    Pres,Prof,Qc=p.read(var,read_adjusted=adj)
	    ii = Pres<=300 ;
# Deve avere almeno 5 records:
	    if len(Prof[ii])>5 :
		NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
		plotmat[:,ip]=NewProf_5m
# FILTRAGGIO per CHLA e N3n dato dai valori troppo alti registrato in superficie:
                if (var_mod=="P_l"):
                   if (plotmat[0,ip] > 0.45):
                      plotmat[:,ip] = np.nan
                if (var_mod=="N3n"): 
                   if (plotmat[0,ip] > 2):
                      plotmat[:,ip] = np.nan
	    timelabel_list.append(p.time)

	    # PLOT FOR THE MODEL
	    TM=MM.modeltime(p)
	    FILENAME = BASEDIR + TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
	    M = readModelProfile(FILENAME,var_mod,p.ID())
	    M_newDepth=np.interp(NewPres_5m,TheMask.zlevels[:max_depth+1],M[:max_depth+1])
	    plotmat_model[:,ip] = M_newDepth


	print var_mod + " " + np.str(len(timelabel_list)) +  p.available_params
        if p.available_params.find(var)<0 : continue
	for ip, p in enumerate(list_float_track):
	    break # esco
	    Lon[ip] = p.lon
	    Lat[ip] = p.lat
	    F=p._my_float
	    Pres,Prof,Profile_adj,Qc = F.read_very_raw(var)
	    Profile_adj
	    ii300 = Pres<=300 ;
	    if len(Prof[ii300])>0 :
		NewProf_5m = np.interp(NewPres_5m,Pres[ii300],Prof[ii300])
	# ESEGUI UN FILTRO GAUSSANO 5M SOPRA E SOTTO IL PUNTO DEFINITO:
		Prof_smooth = meanObj11.compute(NewProf_5m,NewPres_5m)
	    plotmat[:,ip]=NewProf_5m
	#    plotmat[:,ip]=Prof_smooth
	#    xs,ys = np.meshgrid(xlabel_list, dlabels)
	    dt=mpldates.date2num(F.time)
	    timelabel_list.append(dt)

	fig = plt.figure()

	ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2)
	ax2 = plt.subplot2grid((3, 3), (0, 2))
	ax3 = plt.subplot2grid((3, 3), (1, 0), colspan=3)
	ax4 = plt.subplot2grid((3, 3), (2, 0), colspan=3)

	xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_m=ma.masked_invalid(plotmat)
	if (var_mod == 'P_l'):
                quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=0.00,vmax=0.40,cmap="viridis")# default is 'flat'
	if (var_mod == 'O2o'):
		quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=200,vmax=250) #,cmap="jet")# default is 'flat'
        if (var_mod == 'N3n'):
                print len(xs)
                print len(ys)
                print plotmat_m.shape
		quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=0.00,vmax=4) #,cmap="jet")# default is 'flat'
	#Inform matplotlib that the x axis is made by dates
	ax3.set_xlim([T_start2num,T_end2num])
#	ax3.xaxis_date()
        ax3.set_xticklabels([])
	ax3.invert_yaxis()
#	ax3.set_title("FLOAT " + p.name() + ": " + plotvarname[ivar], color = 'b')
	ax3.set_ylabel("depth $[m]$",color = 'k', fontsize=font_s)
	fig.set_size_inches(15,10)
#	fig.set_dpi(300)

	###ax = axs[2]
	xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_model_m=ma.masked_invalid(plotmat_model)
        if (var_mod == 'P_l'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_model_m,shading='flat',vmin=0.00,vmax=0.40,cmap="viridis")# default is 'flat'
        if (var_mod == 'O2o'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_model_m,shading='flat',vmin=200,vmax=250) #,cmap="jet")# default is 'flat'
        if (var_mod == 'N3n'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_model_m,shading='flat',vmin=0.00,vmax=4) #,cmap="jet")# default is 'flat'
	ax4.set_xlim([T_start2num,T_end2num])
	ax4.xaxis_date()
	ax4.invert_yaxis()
#	ax4.set_title("BFM model: " + plotvarname[ivar], color = 'b')
	ax4.set_ylabel("depth $[m]$",color = 'k',fontsize=font_s)
        ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax4.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m")) #("%b%y"))
        for tick in ax4.get_xticklabels():
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
	ind_max_sup=plotmat[0,:].argmax()
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

	fig.savefig(''.join([OUTDIR,'Hov_Float+TRANS_',p.name(),'_',var_mod,'.png']))
	pl.close(fig)
