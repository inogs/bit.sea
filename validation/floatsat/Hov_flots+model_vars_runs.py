import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')


    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--run', '-r',
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
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
import scipy.io.netcdf as NC
import numpy as np
from commons.utils import addsep
import matplotlib.pyplot as pl
from basins import OGS
from validation.online.profileplotter import figure_generator, ncwriter #, add_metadata
import matplotlib
from mhelpers.pgmean import PLGaussianMean
#from profiler_floatsat import *
from instruments.matchup_manager import Matchup_Manager


from commons.mask import Mask
from timeseries.plot import *
meanObj11 = PLGaussianMean(5,1.0)
import matplotlib.pyplot as plt

RUN = args.run

INPUTDIR = '/pico/scratch/userexternal/ateruzzi/' + \
           RUN + '/wrkdir/MODEL/AVE_FREQ_1/'

# output directory, where aveScan.py will be run.

BASEDIR = '/pico/scratch/userexternal/ateruzzi/bit.sea/validation/floatsat/' + \
          RUN + '/PROFILATORE/'


DATESTART = '20150101'
DATE__END = '20150202'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc",filtervar="N1p")

ALL_PROFILES = bio_float.FloatSelector(None,T_INT, Rectangle(-6,36,30,46))


vardescriptorfile="VarDescriptorB.xml"




TheMask=Mask(args.maskfile)
OUTDIR = addsep(args.outdir)

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
VARLIST = ['P_l']
Adj = [True,False,True]
nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11,1.0)

Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

for j in range(0,len(wmo_list)):
    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo_list[j])
    nP = len(list_float_track)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)
    NewPres_5m=np.linspace(0,300,61)
    depths=NewPres_5m
    plotmat = np.zeros([len(depths), len(list_float_track)])
    plotmat_model = np.zeros([len(depths), len(list_float_track)])

    for ivar, var_mod in enumerate(VARLIST):
        timelabel_list = list()
        var = LOVFLOATVARS[var_mod]
        adj=Adj[ivar]

        for ip, p in enumerate(list_float_track):
			Lon[ip] = p.lon
			Lat[ip] = p.lat
			Pres,Prof,Qc=p.read(var,read_adjusted=adj)
			ii = Pres<=300 ;
			if len(Prof[ii])>0 :
				NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
				plotmat[:,ip]=NewProf_5m
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
        if (var_mod == 'P_l'):
		    quadmesh = ax3.pcolormesh(xs, ys, plotmat,shading='flat',vmin=0.00,vmax=0.40)# default is 'flat'
        if (var_mod == 'O2o'):
		    quadmesh = ax3.pcolormesh(xs, ys, plotmat,shading='flat',vmin=200,vmax=250)# default is 'flat'
        if (var_mod == 'N3n'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat,shading='flat',vmin=0.00,vmax=4)# default is 'flat'
		#Inform matplotlib that the x axis is made by dates
        ax3.set_xlim([T_start2num,T_end2num])
        ax3.xaxis_date()
        ax3.invert_yaxis()
        ax3.set_title("FLOAT " + p.name() + ": " + plotvarname[ivar], color = 'b')
        ax3.set_ylabel("depth $[m]$",color = 'k')
        fig.set_size_inches(15,10)
        fig.set_dpi(300)

		###ax = axs[2]
        xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        if (var_mod == 'P_l'):
		    quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=0.00,vmax=0.40)# default is 'flat'
        if (var_mod == 'O2o'):
		    quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=200,vmax=250)# default is 'flat'
        if (var_mod == 'N3n'):
		    quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=0.00,vmax=4)# default is 'flat'
        ax4.set_xlim([T_start2num,T_end2num])
        ax4.xaxis_date()
        ax4.invert_yaxis()
        ax4.set_title(RUN + " " + plotvarname[ivar], color = 'b')
        ax4.set_ylabel("depth $[m]$",color = 'k')

        colorbar_bottom= ax4.get_position().ymin
        colorbar_extent= ax3.get_position().ymax - colorbar_bottom

	#      cbaxes = fig.add_axes([0.95, 0.01, 0.02, 0.95]) 
        cbaxes = fig.add_axes([0.93, colorbar_bottom , 0.02, colorbar_extent])
        fig.colorbar(quadmesh, cax = cbaxes, )  


        c_lon,c_lat=coastline.get()

        ax1.plot(c_lon,c_lat,'k')
        ax1.plot(Lon,Lat,'r.')
        ax1.set_title("TRAJECTORY of FLOAT " + p.name() + " - " + var_mod, color = 'r')
        ind_max_sup=plotmat[0,:].argmax()
		#      print Lon[ind_max_sup],Lat[ind_max_sup]
		#	ax1.plot(Lon[ind_max_sup],Lat[ind_max_sup],'g.')
        ax1.plot(Lon[0],Lat[0],'bo')
        ax1.plot(Lon[0],Lat[0],'bx')
        ax1.set_xlim([-10,36])
        ax1.set_ylabel("LAT",color = 'k')
        ax1.set_xlabel("LON",color = 'k')

        extent=4 #degrees
        ipp = len(list_float_track)
        ax2.plot(c_lon,c_lat,'k')
        ax2.plot(Lon,Lat,'ro')
        #	ax2.plot(Lon[ind_max_sup],Lat[ind_max_sup],'go')
        ax2.plot(Lon[0],Lat[0],'bo')
        ax2.plot(Lon[0],Lat[0],'bx',markersize=12)
        ax2.set_xlim([np.min(Lon[:ipp]) -extent/2, np.max(Lon[:ipp]) +extent/2])
        ax2.set_ylim([np.min(Lat[:ipp]) -extent/2, np.max(Lat[:ipp]) +extent/2])

        print('-----------------------')

        fig.savefig(''.join([OUTDIR,'Hov_Float+TRANS_',p.name(),'_',var_mod,'.png']))
        pl.close(fig)
