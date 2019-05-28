import argparse
from commons.utils import addsep
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')


    parser.add_argument(   '--indir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)

import os,sys
import scipy.io.netcdf as NC
from commons.mask import Mask
from instruments import matchup_manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Region, Rectangle
from layer_integral import coastline
import instruments
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
import scipy.io.netcdf as NC
import numpy as np
from commons.utils import addsep
import pylab as pl
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

from commons.layer import Layer
from scipy.stats import pearsonr
TheMask=Mask(args.maskfile)

font_s =  15 
font_s2 = 3
label_s = 15

layer=Layer(0,300)

#######
# For second plot:

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from commons.utils import addsep
from profiler import ALL_PROFILES,TL,BASEDIR
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
import basins.V2 as OGS
from datetime import datetime
from datetime import timedelta

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','dNit_dz']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
MED_PROFILES = bio_float.FloatSelector(None,TL.timeinterval,OGS.med)
wmo_list=bio_float.get_wmo_list(MED_PROFILES)

izmax = TheMask.getDepthIndex(200) 

################

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

def get_level_depth(TheMask,lev):

    i0 = TheMask.getDepthIndex(lev)
    i1 = i0 + 1
      
    diff0 = lev - TheMask.zlevels[i0]
    diff1 = lev - TheMask.zlevels[i1]
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

max_depth = get_level_depth(TheMask,300)

MM = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
plotvarname = [r'Chl $[mg/m^3]$',r'Oxy $[mmol/m^3]$',r'Nitr $[mmol/m^3]$'] #,r'Temp $[^\circ C]$','Sal']
VARLIST = ['P_l','O2o','N3n']
Adj = [True,False,True]
nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11,1.0)

Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

print wmo_list

for j in range(0,len(wmo_list)):
      list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo_list[j])
      nP = len(list_float_track)
      Lon = np.zeros((nP,), np.float64)
      Lat = np.zeros((nP,), np.float64)

      for ivar, var_mod in enumerate(VARLIST):
        max_depth = get_level_depth(TheMask,300)
        NewPres_5m=np.linspace(0,300,121)
        bt=300

        depths=NewPres_5m
        plotmat = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_model = np.zeros([len(depths), len(list_float_track)])*np.nan
	timelabel_list = list()
        var = LOVFLOATVARS[var_mod]
        adj=Adj[ivar]
        
	for ip, p in enumerate(list_float_track):
	    Lon[ip] = p.lon
	    Lat[ip] = p.lat
	    Pres,Prof,Qc=p.read(var,read_adjusted=adj)
	    ii = Pres<=bt ;
            if len(Prof[ii])>1 : # At least two records
		NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
		plotmat[:,ip]=NewProf_5m
                if (var_mod=="P_l"):
                    # FILTER OUT PROFILES WITH HIGH VALUES AT SURF
                    if (plotmat[0,ip] > 0.8):
                      plotmat[:,ip] = np.nan
#                      plotmat_model[:,ip] = np.nan
                if (var_mod=="N3n"): 
                   if (plotmat[0,ip] > 2):
                      plotmat[:,ip] = np.nan
#                      plotmat_model[:,ip] = np.nan
	    timelabel_list.append(p.time)

	    # PLOT FOR THE MODEL
	    TM=MM.modeltime(p)
	    FILENAME = BASEDIR + TM.strftime("PROFILES/ave.%Y%m%d-12:00:00.profiles.nc")
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
	    ii300 = Pres<=bt;
	    if len(Prof[ii300])>0 :
		NewProf_5m = np.interp(NewPres_5m,Pres[ii300],Prof[ii300])
	# ESEGUI UN FILTRO GAUSSANO 5M SOPRA E SOTTO IL PUNTO DEFINITO:
		Prof_smooth = meanObj11.compute(NewProf_5m,NewPres_5m)
	    plotmat[:,ip]=NewProf_5m
	    dt=mpldates.date2num(F.time)
	    timelabel_list.append(dt)

	fig = plt.figure()

        # AXES FOR MAPS
	ax1 = plt.subplot2grid((7, 3), (0, 0), colspan=2)
	ax2 = plt.subplot2grid((7, 3), (0, 2))
        # AXES PLOT THE HOVMOELLER
	ax3 = plt.subplot2grid((7, 3), (1, 0), colspan=3)
	ax4 = plt.subplot2grid((7, 3), (2, 0), colspan=3)
        # AXES FOR STATISTICS
        ax5 = plt.subplot2grid((7, 3), (3, 0), colspan=3)
        ax6 = plt.subplot2grid((7, 3), (4, 0), colspan=3)
        ax7 = plt.subplot2grid((7, 3), (5, 0), colspan=3)
        ax8 = plt.subplot2grid((7, 3), (6, 0), colspan=3)

	xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_m=ma.masked_invalid(plotmat)
        # PLOT HOVMOELLER OF FLOAT
	if (var_mod == 'P_l'):
                quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=0.00,vmax=0.40,cmap="viridis")# default is 'flat'
	if (var_mod == 'O2o'):
                # SET RANGES BETWEEN 160 AND 250
                quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=160,vmax=250) #,cmap="jet")# default is 'flat'
        if (var_mod == 'N3n'):
                print len(xs)
                print len(ys)
                print plotmat_m.shape
		quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=0.00,vmax=4) #,cmap="jet")# default is 'flat'
	#Inform matplotlib that the x axis is made by dates
        ax3.set_xlim([T_start2num,timelabel_list[-1]])  # SET THE END BY LENGHT OF TIMESERIES
#	ax3.xaxis_date()
#        ax3.set_xticklabels([])
#	ax3.invert_yaxis()
#	ax3.set_title("FLOAT " + p.name() + ": " + plotvarname[ivar], color = 'b')
	ax3.set_ylabel("depth $[m]$",color = 'k', fontsize=font_s)
        fig.set_size_inches(12,15)
	fig.set_dpi(300)

	###ax = axs[2]
        # PLOT HOVMOELLER OF MODEL
	xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_model_m=ma.masked_invalid(plotmat_model)
        if (var_mod == 'P_l'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=0.00,vmax=0.40,cmap="viridis")# default is 'flat'
        if (var_mod == 'O2o'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=160,vmax=250) #,cmap="jet")# default is 'flat'
        if (var_mod == 'N3n'):
                quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=0.00,vmax=4) #,cmap="jet")# default is 'flat'
#	ax4.set_xlim([T_start2num,T_end2num])
        ax3.set_xlim([T_start2num,timelabel_list[-1]])
        ax3.invert_yaxis()
        ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))


        ax4.set_xlim([T_start2num,timelabel_list[-1]])
#	ax4.xaxis_date()

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
###	cbar=fig.colorbar(quadmesh, cax = cbaxes, )  
###	ticklabs = cbar.ax.get_yticklabels()
###	cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)

        # MAP OF TRAJECTORY:
	c_lon,c_lat=coastline.get()

	ax1.plot(c_lon,c_lat,'k')
	ax1.plot(Lon,Lat,'r.',markersize=font_s2)
#	ax1.set_title("TRAJECTORY of FLOAT " + p.name() + " - " + var_mod, color = 'r', fontsize=font_s)
        ax1.set_title("TRAJECTORY of FLOAT " + p.name() , color = 'r', fontsize=font_s)
	ind_max_sup=plotmat[0,:].argmax()
        ax1.plot(Lon[0],Lat[0],'bo',markersize=font_s2)
	ax1.plot(Lon[0],Lat[0],'bx')
	ax1.set_xlim([-5.5,36])
#	ax1.set_ylabel("LAT",color = 'k',fontsize=font_s)
#	ax1.set_xlabel("LON",color = 'k')
	extent=4 #degrees
	ipp = len(list_float_track)
	ax2.plot(c_lon,c_lat,'k')
	ax2.plot(Lon,Lat,'ro',markersize=font_s2)
#	ax2.plot(Lon[ind_max_sup],Lat[ind_max_sup],'go')
	ax2.plot(Lon[0],Lat[0],'bo',markersize=font_s2)
	ax2.plot(Lon[0],Lat[0],'bx',markersize=font_s)
	ax2.set_xlim([np.min(Lon[:ipp]) -extent/2, np.max(Lon[:ipp]) +extent/2])
	ax2.set_ylim([np.min(Lat[:ipp]) -extent/2, np.max(Lat[:ipp]) +extent/2])



# SECOND FIGURE ON STATISTICS:

        izmax = TheMask.getDepthIndex(200) # Max Index for depth 200m

#        INPUT_FILE = INDIR + wmo + ".nc"
        # READ STATISTICS FILE *.nc
        INPUT_FILE = INDIR + wmo_list[j] + ".nc"
        print INPUT_FILE
        A = ncreader(INPUT_FILE)
        times = timelabel_list
        print var_mod
#        OUTFILE = OUTDIR + var + "_" + wmo + ".png"
        OUTFILE = OUTDIR + var_mod + "_" + wmo_list[j] + ".png"
        print OUTFILE

	model, ref =A.plotdata(var_mod,'Int_0-200')
        surf_model, surf_ref = A.plotdata(var_mod,'SurfVal')
#	if (~np.isnan(model).all() == True) or (~np.isnan(ref).all() == True): 
        if (1 == 1):
            ax6.plot(times,  ref,'.b',label='REF INTEG')
            ax6.plot(times,model,'b',label='MOD INTEG')
	    
	    ax5.plot(times,  surf_ref,'.b',label='REF SURF')
            ax5.plot(times,  surf_model,'b',label='MOD SURF')
	    if ( var_mod == "P_l" ):
	        ax6.set_ylabel('INTG 0-200 \n $[mg{\  } m^{-3}]$',fontsize=15)
                ax5.set_ylabel('SURF \n $[mg{\  } m^{-3}]$',fontsize=15)
                ax5.set_ylim(0,0.5)
                ax6.set_ylim(0,0.22)
#                ax7.set_ylim(40,200)
#                ax7.set_ylim(0,200)
#                ax8.set_ylim(0,200)
	    if ( var_mod == "O2o" ):
                ax6.set_ylabel('INTG \n $[mmol{\  } m^{-3}]$',fontsize=15)
                ax5.set_ylabel('SURF \n $[mmol{\  } m^{-3}]$',fontsize=15)
            if ( var_mod == "N3n" ):
                ax6.set_ylabel('INTG 0-350m \n $[mmol{\  } m^{-3}]$',fontsize=15)
#		axes[2].set_ylabel('INTEG 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
#                axes[2].set_ylim(0,8)
                ax5.set_ylim(0,1.8)
                ax6.set_ylim(0,2.8)
#                axes[4].set_ylim(40,200)
                ax5.set_ylabel('SURF \n $[mmol{\  } m^{-3}]$',fontsize=15)
                ax8.set_ylim(0,350) # RANGE FOR THE NITRICLINO
#                ax5.set_ylabel('at 50m depth \n $[mmol{\  } m^{-3}]$',fontsize=15)
#	    legend = axes[2].legend(loc='upper left', shadow=True, fontsize=12)
 	    model_corr , ref_corr =A.plotdata(var_mod,'Corr')
	    times_r = times
	    for icr , cr in enumerate(ref_corr): 
                print cr
                if ( ( ~np.isnan(cr) <= 0 ) | np.isnan(cr) ):
		    ref_corr[icr] = np.nan
#                    times_r.remove(times[icr])

#	    ax7.plot(times,ref_corr,'b')
#	    ax7.set_ylabel('CORR',fontsize=15)

	    ax2.set_xticklabels([])
	    ax3.set_xticklabels([])
            ax4.set_xticklabels([])
            ax5.set_xticklabels([])
            ax6.set_xticklabels([])
            if (np.isnan(ref_corr).all() == False ):
                ax7.plot(times,ref_corr,'b')
                ax7.set_ylabel('CORR',fontsize=15)
                ax7.set_xticklabels([])
                ax7.set_ylim(0,1)

                ax7.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
                ax7.tick_params(axis = 'both', which = 'major', labelsize = label_s)
                ax7.set_xlim([T_start2num,timelabel_list[-1]])


        if (var_mod == "P_l"): 
            model_dcm, ref_dcm =A.plotdata(var_mod,'DCM')
	    model_mld, ref_mld =A.plotdata(var_mod,'z_01')
#            if (model_mld > 150):
#	       model_mld = np.nan
 
#	    if (~np.isnan(model_dcm).all() == True) or (~np.isnan(ref_dcm).all() == True):
            if (1 == 1):
                ax8.invert_yaxis()
                ax8.plot(times,  ref_dcm,'.b',label='DCM REF')
                ax8.plot(times,model_dcm,'b',label='DCM MOD')
                ax8.plot(times, ref_mld,'.r',label='MWB REF')
                ax8.plot(times,model_mld,'r',label='MWB MOD')
#		axes[4].plot(times,np.ones_like(times)* np.nanmean(ref_dcm),'r',linewidth=3) 
#		axes[4].plot(times,np.ones_like(times)* np.nanmean(model_dcm),'b',linewidth=3) #marker='.')

	        ax8.set_ylabel('DCM/MWB $[m]$',fontsize=15)
#		axes[4].xaxis_date()
#		axes[4].xaxis.set_major_locator(mdates.MonthLocator())
		ax8.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
                ax8.set_xlim([T_start2num,timelabel_list[-1]])
                ax8.set_ylim([200,0])


		xlabels = ax8.get_xticklabels()
		plt.setp(xlabels, rotation=30,fontsize=15)


        if (var_mod == "N3n"):
            model_nit, ref_nit =A.plotdata(var_mod,'Nit_1')
            print len(model_nit)
            model_nit2, ref_nit2 = A.plotdata(var_mod,'dNit_dz')
            print [ ref_nit2[ ~np.isnan(ref_nit2) ] <= 10.0 ]
            ref_nit2[ np.isnan(ref_nit2) ] = 0.0
            jj=[ref_nit2[ ~np.isnan(ref_nit2) ] <= 10.0 ] 
            ref_nit2[jj] = np.nan
            print ref_nit2
	    if (~np.isnan(ref_nit).all() == True) or (~np.isnan(model_nit).all() == True):
            	ax8.plot(times,  ref_nit,'.b',label='REF')
            	ax8.plot(times,model_nit,'b',label='MOD')
            	ax8.invert_yaxis()
	    	ax8.set_ylabel('NITRCL 1/2 $[m]$',fontsize=15)
#                pint ref_nit2
                ax8.plot(times, ref_nit2,'.r',label='dNit REF')
                ax8.plot(times,model_nit2,'r',label='dNit MOD') 
            else: 
		ax8.plot(times,  np.ones_like(times))
            ax8.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = ax8.get_xticklabels()
            plt.setp(xlabels, rotation=30)
#	    legend = axes[4].legend(loc='upper left', shadow=True, fontsize=12)

        if (var_mod == "O2o"):
	    ax8.plot(times,  np.ones_like(times))
            ax8.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            xlabels = ax8.get_xticklabels()
            plt.setp(xlabels, rotation=30,fontsize=15)

        cbar=fig.colorbar(quadmesh, cax = cbaxes, )
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)

        ax5.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax6.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
#        ax7.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))  # moved above
        ax8.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))

        ax3.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax3.set_xlim([T_start2num,timelabel_list[-1]])
        ax4.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax4.set_xlim([T_start2num,timelabel_list[-1]])

        ax5.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax5.set_xlim([T_start2num,timelabel_list[-1]])
        ax6.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax6.set_xlim([T_start2num,timelabel_list[-1]])
#        ax7.tick_params(axis = 'both', which = 'major', labelsize = label_s)   # moved above
#        ax7.set_xlim([T_start2num,timelabel_list[-1]])                         # moved above
        ax8.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax8.set_xlim([T_start2num,timelabel_list[-1]])
        ax5.xaxis.grid(True)
        ax6.xaxis.grid(True)
        ax7.xaxis.grid(True)
        ax8.xaxis.grid(True)


        fig.savefig(OUTFILE)
	plt.close(fig)
