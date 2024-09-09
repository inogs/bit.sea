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

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help = "Directory where metrics are saved")

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
import scipy.io.netcdf as NC

from commons.layer import Layer
from commons.mask import Mask
from commons.utils import addsep
from metrics import *

from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from instruments.matchup_manager import Matchup_Manager

from profilerPhys_comparison2017 import *

from layer_integral import coastline
from basins import OGS

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from timeseries.plot import *
import numpy.ma as ma
from mhelpers.pgmean import PLGaussianMean

from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader

TheMask = Mask(args.maskfile)

OUTDIR = addsep(args.outdir)

DIRSTATS = {}
for run in [RUN_REF,RUN_DA]:
	DIRSTATS[run] = addsep(args.indir) + '/' + run + '/STATS/'


font_s =  15 
label_s = 15

def my_Hovmoeller_diagram(plotmat, xs,ys, fig=None, ax=None):
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
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

max_depth = get_level300(TheMask)
layerStats = Layer(0,200)

# plotvarname = [r'Chl $[mg/m^3]$',r'Oxy $[mmol/m^3]$',r'Nitr $[mmol/m^3]$'] #,r'Temp $[^\circ C]$','Sal']
VARLIST = ['votemper','vosaline']
Adj = {
	'votemper':  False,
	'vosaline': True,
	'O2o':  False,
	'N3n':  True,
    }
nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11,1.0)

Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

MM = {}
for run in [RUN_REF,RUN_DA]:
    MM[run] = Matchup_Manager(ALL_PROFILES,TL,BASEDIR[run])


#print wmo_list

DICTstyle = {
    RUN_REF: ['-k','-g','-b'],
    RUN_DA:  ['--k','--g','--b'],
    'Float': [':k',':g',':b'],
}

for j,wmo in enumerate(wmo_list):
# for j,wmo in enumerate([wmo_list[0]]):
    print(wmo)
    stats = {}
    for run in [RUN_REF,RUN_DA]:
		filestats = DIRSTATS[run] + '/' + wmo + 'phys.nc'
		stats[run] = ncreader(filestats)

    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    nP = len(list_float_track)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)
    NewPres_5m=np.linspace(0,300,61)
    depths=NewPres_5m

    for ivar, var_mod in enumerate(VARLIST):
#     if (ivar==0):
        plotmat = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_model = {}
        for run in [RUN_REF,RUN_DA]:
			plotmat_model[run] = np.zeros([len(depths), len(list_float_track)])*np.nan
        timelabel_list = list()
        var = LOVFLOATVARS[var_mod]
        adj=Adj[var_mod]

        for ip, p in enumerate(list_float_track):
            Lon[ip] = p.lon
            Lat[ip] = p.lat
            Pres,Prof,Qc = p.read(var,read_adjusted=adj)
            ii = Pres<=300
			# Deve avere almeno 5 records:
            if len(Prof[ii])>5 :
				NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
				plotmat[:,ip]=NewProf_5m
				# FILTRAGGIO per CHLA e N3n dato dai valori troppo alti registrato in superficie:
				# if (var_mod=="P_l"):
				# 	if (plotmat[0,ip] > 0.45):
				# 		plotmat[:,ip] = np.nan
				# if (var_mod=="N3n"): 
				# 	if (plotmat[0,ip] > 2):
				# 		plotmat[:,ip] = np.nan
            timelabel_list.append(p.time)

# PLOT FOR THE MODEL
            for run in [RUN_REF,RUN_DA]:
                TM = MM[run].modeltime(p)
                FILENAME = BASEDIR[run] + TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
                M = readModelProfile(FILENAME,var_mod,p.ID())
                M_newDepth=np.interp(NewPres_5m,TheMask.zlevels[:max_depth+1],M[:max_depth+1])
                plotmat_model[run][:,ip] = M_newDepth


        print var_mod + " " + np.str(len(timelabel_list)) +  p.available_params

    	plt.close('all')
        fig,axs = plt.subplots(2,1,sharex=True,figsize=[15,10])

        ax1 = axs[0]
        plt.sca(ax1)
        plt.title(wmo + ' RunDAphys - RunHCphys')

        xs,ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_diffmod_m = ma.masked_invalid(plotmat_model[RUN_DA]-plotmat_model[RUN_REF])
        if (var_mod == 'votemper'):
            quadmesh = ax1.pcolormesh(xs, ys, plotmat_diffmod_m, \
                shading='flat',vmin=-1,vmax=1,cmap="coolwarm")
        if (var_mod == 'vosaline'):
			quadmesh = ax1.pcolormesh(xs, ys, plotmat_diffmod_m, \
				shading='flat',vmin=-.5,vmax=.5,cmap="coolwarm")

        ax1.set_xlim([T_start2num,T_end2num])
        ax1.xaxis_date()
        ax1.invert_yaxis()
		#ax1.set_title("BFM model: " + plotvarname[ivar], color = 'b')
        ax1.set_ylabel("depth $[m]$",color = 'k',fontsize=font_s)
        ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m")) #("%b%y"))
        for tick in ax1.get_xticklabels():
            tick.set_rotation(45)

        ax1.tick_params(axis = 'both', which = 'major', labelsize = label_s)

        colorbar_bottom= ax1.get_position().ymin
        colorbar_extent= ax1.get_position().ymax - colorbar_bottom
        cbaxes = fig.add_axes([0.93, colorbar_bottom , 0.02, colorbar_extent])
        cbar = fig.colorbar(quadmesh, cax = cbaxes, )
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)



        ax2 = axs[1]

        if (var_mod=='votemper'):
            for run in [RUN_REF,RUN_DA]:
                model_tcl1,ref_tcl1 = stats[run].plotdata(var_mod,'Termocl1')
                model_tcl2,ref_tcl2 = stats[run].plotdata(var_mod,'Termocl2')
                ax2.plot(mpldates.date2num(timelabel_list),model_tcl1, \
                        DICTstyle[run][0],label='14 C ' + run)
                ax2.plot(mpldates.date2num(timelabel_list),model_tcl2, \
                        DICTstyle[run][1],label='16 C ' + run)
            ax2.plot(mpldates.date2num(timelabel_list),ref_tcl1, \
                    DICTstyle['Float'][0],label='Float 14 C')
            ax2.plot(mpldates.date2num(timelabel_list),ref_tcl2, \
                    DICTstyle['Float'][1],label='Float 16 C')
            ax2.invert_yaxis()
            ax2.set_ylabel("Depth $[m]$",color = 'k',fontsize=font_s)
            ax2.tick_params(axis='both', which='major',labelsize=label_s)
            
            plt.sca(ax2)
            plt.grid()
            plt.legend(loc='best')
        

        fig.savefig(''.join([OUTDIR,'HovdiffPhys_',p.name(),'_',var_mod,'.png']))


plt.show(block=False)
