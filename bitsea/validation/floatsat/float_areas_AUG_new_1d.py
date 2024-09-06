import os,sys
import scipy.io.netcdf as NC
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
import scipy.io.netcdf as NC
import numpy as np
import matplotlib.pyplot as pl
from commons.mask import Mask
from postproc import maskload   #read_Positions_for_Pointprofiles
from instruments import matchup_manager #	readModelProfile
from datetime import datetime
import commons.genUserDateList as DL


from basins.region import Rectangle
import netCDF4 as nc
import csv
import matplotlib.dates as mpldates

pl.close('all')

SCRATCHDIR = '/pico/scratch/userexternal/ateruzzi/'
RUN = 'SAT_01'
maskfile = SCRATCHDIR + '/DA_FLOAT_SAT/Summer/RUN_' + \
             RUN +'/wrkdir/MODEL/meshmask.nc'
TheMask=Mask(maskfile)
PUNTIFILE='punti_20150804_20150811_20150818_20150825.dat'


# ASSIM at 12
#INPUTDIR_model='/pico/scratch/userexternal/lmariott/FLOAT_DA_01/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/'

#start = datetime.datetime(2015,1,6)
#end   = datetime.datetime(2015,2,6)
#step = datetime.timedelta(days=5)

#ass_dates_F = []


# ASSIM at 00
#INPUTDIR_model='/pico/scratch/userexternal/lmariott/FLOAT_DA_03_newstd/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/'
INPUTDIR_model = SCRATCHDIR + 'DA_FLOAT_SAT/Summer/RUN_' + RUN + \
                 '/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/'
INPUTDIR_hind = SCRATCHDIR + 'DA_FLOAT_SAT/Summer/RUN_CR' + \
                 '/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/'
#INPUTDIR='/pico/home/userexternal/lfeudale/work/DA/test1/AUG/PROFILES/'
#OUTDIR  = '/pico/scratch/userexternal/lfeudale/DA_validation/output/FLOAT_DA_03_newstd/'
OUTDIR = 'DA_FLOAT_SAT/Summer/RUN_' + RUN + '/'

def my_Hovmoeller_diagram(plotmat, xs,ys, fig=None, ax=None):
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
    quadmesh = ax.pcolormesh(xs, ys, plotmat,shading='gouraud')# default is 'flat'
    #Inform matplotlib that the x axis is made by dates
    ax.xaxis_date()
    ax.invert_yaxis()
    return fig, ax, quadmesh

def readModelProfile(filename,var, wmo):
        '''
        Reads profiles produced by aveScan.py.
        In these files each variable has dimensions (jpk, nProfiles)
        And each profile is identified by the corresponding wmo
        '''
        ncIN = NC.netcdf_file(filename,'r')

        M = ncIN.variables[var].data.copy()
        iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
        ncIN.close()
        Profile = M[iProfile,:]

        return Profile

#TI1 = TimeInterval('20150801','20150901','%Y%m%d')
TI1 = TimeInterval('20150801','20150909','%Y%m%d')

# ASSIM at 12
#ass_dates = ['20150107','20150114','20150121','20150128']
# ASSIM at 00
ass_dates = ['20150804','20150811','20150818','20150825','20150901']
#ass_dates_F = [datetime(2015,01,06), datetime(2015,01,13), datetime(2015,01,20), datetime(2015,01,27), datetime(2015,02,03)]
#assim_dates=[]
#for j in [0,1,2,3,4]: assim_dates.append(datetime.strptime(ass_dates[j]+"00","%Y%m%d%H"))
#assim_dates_F = TimeList.fromfilenames(TI1,INPUTDIR_model,"RST.*P_l.nc",prefix='RST.')
import commons.genUserDateList as DL

#assim_dates_F =  DL.getTimeList("20150106-00:00:00","20150209-00:00:00","days=5")
assim_dates_F =  DL.getTimeList("20150804-00:00:00","20150907-00:00:00","days=7")


varname = ['CHLA']
plotvarname = [r'Chl $[mg/m^3]$']
model_varname = 'P_l'
#        ref_varname = M.reference_var(p, model_varname)

Points = maskload.read_Positions_for_Pointprofiles(PUNTIFILE)
max_depth = 26
extragrid = 2./16 # -> 2 grid_points = 2 * 1/16

print Points['Name']

#TL = TimeList.fromfilenames(TI1, INPUTDIR,"ave*.nc",filtervar="nc")
TL = TimeList.fromfilenames(TI1, INPUTDIR_model,"ave*.P_l.nc",filtervar="P_l",dateformat='%Y%m%d')
new_TL = [datetime(t.year,t.month,t.day,0)   for t in TL.Timelist  ]
ndays=len(TL.Timelist)

wmoset = set()
for j, jj in enumerate(Points['Name']):
	wmo=Points['Name'][j][:7]
	wmoset.add(wmo)

wmolist=list(wmoset)

table_integ = np.zeros([1 + len(wmolist), 1 + len(TL.Timelist)])
table_integ2= np.zeros([1 + len(wmolist), 1 + len(TL.Timelist)])
#for i in np.arange(0,39): 
for i in np.arange(0,ndays):
	table_integ[0,i+1]=TL.Timelist[i].strftime("%Y%m%d")
	table_integ2[0,i+1]=TL.Timelist[i].strftime("%Y%m%d")

for iwmo, st_wmo in enumerate(wmolist):
        print st_wmo
#	fig , axs = pl.subplots(3,2)
        
        Point_list=[]	
	Point_LAT=[]
	Point_LON=[]
	for ii, s in enumerate(Points['Name']):
   	     if (s.find(st_wmo)==0) : 
                 Point_list.append(s)
		 Point_LAT.append(Points['Lat'][ii])
		 Point_LON.append(Points['Lon'][ii])
	print len(Point_list)

	print "MINLON " + str(min(Point_LON))
        print "MAXLON " + str(max(Point_LON))
	print "MINLAT " + str(min(Point_LAT))
        print "MAXLAT " + str(max(Point_LAT))

	minlon = min(Point_LON) - extragrid
	maxlon = max(Point_LON) + extragrid
	minlat = min(Point_LAT) - extragrid
        maxlat = max(Point_LAT) + extragrid
	
	area_float = Rectangle(minlon,maxlon,minlat,maxlat)
	
	plotmat = np.zeros([max_depth, len(TL.Timelist)])
	plotmat_hind = np.zeros([max_depth, len(TL.Timelist)])
	plotmat_diff = np.zeros([max_depth, len(TL.Timelist)])
        plotmat_fr   = np.zeros([max_depth, len(TL.Timelist)])

	table_integ[iwmo+1,0] = wmolist[iwmo]
	table_integ2[iwmo+1,0] = wmolist[iwmo]
	for ip, p in enumerate(new_TL):  #.Timelist):
#		FILENAME = INPUTDIR_model + p.strftime("ave.%Y%m%d-%H:%M:%S.") + model_varname + ".nc"
                FILENAME = INPUTDIR_model + p.strftime("ave.%Y%m%d-12:00:00.") + model_varname + ".nc"
#		FILENAME_hind = INPUTDIR_hind + p.strftime("ave.%Y%m%d-%H:%M:%S.") + model_varname + ".nc"
                FILENAME_hind = INPUTDIR_hind + p.strftime("ave.%Y%m%d-12:00:00.") + model_varname + ".nc"
#                print FILENAME
		f = nc.Dataset(FILENAME)
		f_hind =  nc.Dataset(FILENAME_hind)
		Chl_field = f.variables["P_l"]
		Chl_field_hind = f_hind.variables["P_l"]
		lon = f.variables["lon"]
                lat = f.variables["lat"]
		
		# select the regions
		latidx = (lat >= minlat) & (lat <=maxlat)
		lonidx = (lon >= minlon) & (lon <=maxlon)

	#	Chl = Chl_field[:, latidx][..., lonidx]
		Chl_mean_profile =  np.mean(np.mean(Chl_field[:,:max_depth, latidx,lonidx],axis=3),axis=2)
		plotmat[:,ip] = Chl_mean_profile

		Chl_mean_profile_hind =  np.mean(np.mean(Chl_field_hind[:,:max_depth, latidx,lonidx],axis=3),axis=2)
                plotmat_hind[:,ip] = Chl_mean_profile_hind

		plotmat_diff[:,ip] = Chl_mean_profile - Chl_mean_profile_hind
        st_d = 3 #number of days before the first assimilation
	tstep = 7 #time step of the assimilation
#	wcount=(np.arange(0,32,1)-7)/7
        wcount=(np.arange(0,ndays,1)-st_d)/tstep
#	for j in np.arange(7,32,1): 
        for j in np.arange(st_d,ndays,1):
#        for j in np.arange(st_d,36,1):
	        	
		persis= abs(plotmat_diff[:22,j] - plotmat_diff[:22,st_d+tstep*wcount[j]]) # w is the week of reference
		persis2 = abs(plotmat_diff[:22,j] / plotmat_diff[:22,st_d+tstep*wcount[j]])
		#print table_integ2[0,j+1]
		#print persis2
#		z_integ = np.sum(persis*np.diff(TheMask.dz[:23])) 
                z_integ = np.sum(persis*np.diff((np.append(0,TheMask.zlevels[:22]))))/TheMask.zlevels[21]
                z_integ2= np.sum(persis2*np.diff((np.append(0,TheMask.zlevels[:22]))))/TheMask.zlevels[21]
		#print z_integ, z_integ2
		table_integ[iwmo+1,j+1] = z_integ
		table_integ2[iwmo+1,j+1] = z_integ2
		import sys
#		sys.exit()
	

#        dt = mpldates.date2num(TL.Timelist)
        dt = mpldates.date2num(new_TL)
        xs,ys = np.meshgrid(dt, TheMask.zlevels[:max_depth])
        fig , axs = pl.subplots(3,1)
        quadmesh = axs[0].pcolormesh(xs, ys, plotmat,shading='flat',vmin=0.0,vmax=0.3)
	fig.colorbar(quadmesh,ax=axs[0])
	axs[0].xaxis_date()
#	xfmt = ass_dates_F.DateFormatter('%d-%m-%y')
#        ax.xaxis.set_major_formatter(xfmt)
#	axs[2].xaxis.set_ticks(ass_dates_F)
#	axs[2].set_xticklabels(ass_dates_F, rotation=20, horizontalalignment = 'right')
        quadmesh_hind = axs[1].pcolormesh(xs, ys, plotmat_hind,shading='flat',vmin=0.0,vmax=0.3)
	fig.colorbar(quadmesh_hind,ax=axs[1])

        quadmesh_diff = axs[2].pcolormesh(xs, ys, plotmat_diff,shading='flat',vmin=-0.03,vmax=0.03, cmap=pl.cm.get_cmap('seismic'))
	fig.colorbar(quadmesh_diff)

	title = st_wmo + "_AUG_" + RUN + varname[0] + "_" + str(min(Point_LON)) + "_" + str(max(Point_LON)) + "_"  + str(min(Point_LAT)) + "_" + str(max(Point_LAT))
#	axs[0].set_title( st_wmo + " - " + plotvarname[0] + "[" + str(min(Point_LON)) + "-" + str(max(Point_LON)) + ","  + str(min(Point_LAT)) + "-" + str(max(Point_LAT)) + "]")

        axs[0].text(axs[0].get_xlim()[1]+3,axs[0].get_ylim()[0]-3,r'$\mathrm{ [\, mg \, CHLa /m^3]} $', \
	        fontsize=12,horizontalalignment='right',verticalalignment='bottom')
	axs[0].text(axs[0].get_xlim()[1]-3,axs[0].get_ylim()[0],"RUN_" + RUN,\
	        fontsize=14,horizontalalignment='right',verticalalignment='bottom')
        axs[1].text(axs[1].get_xlim()[1]-3,axs[1].get_ylim()[0],"CR", \
	        fontsize=14,horizontalalignment='right',verticalalignment='bottom')
        axs[2].text(axs[2].get_xlim()[1]-3,axs[2].get_ylim()[0],"RUN_" + RUN + "-CR", \
	        fontsize=14,horizontalalignment='right',verticalalignment='bottom')

#	axs[1].set_title("Hindcast")
#        axs[2].set_title("Assimilation-Hindcast")

	for k in [0,1,2]:
        	axs[k].xaxis_date()
        	axs[k].invert_yaxis()
		axs[k].set_ylabel('depth $[m]$')
		axs[k].set_xlim(axs[k].set_xlim(dt[0],dt[-1]))
		axs[k].set_ylim(axs[k].set_ylim(TheMask.zlevels[25],TheMask.zlevels[0]))
# TheMask.zlevels[25] --> 209.48465m
# TheMask.zlevels[21] --> 148.32494
		axs[k].tick_params(axis='both', labelsize=10)
		axs[k].grid(axis='both')
#		if (k != 1):
#        	    for j in [0,1,2,3,4]: 
# 			axs[k].plot([mpldates.date2num(assim_dates_F[j]),mpldates.date2num(assim_dates_F[j])],[TheMask.zlevels[25],TheMask.zlevels[0]],'k--', lw=2)
                for j in np.arange(0,len(assim_dates_F)): 
		    axs[k].plot([mpldates.date2num(assim_dates_F[j]), \
		                 mpldates.date2num(assim_dates_F[j])], \
				 [TheMask.zlevels[25],TheMask.zlevels[0]], \
				 'k--', lw=2)
        axs[0].set_xticklabels([])
        axs[1].set_xticklabels([])
        axs[2].xaxis.set_ticks(assim_dates_F)
        axs[2].xaxis.set_major_formatter(mpldates.DateFormatter("%d-%m-%Y"))
        xlabels = axs[2].get_xticklabels()
        pl.setp(xlabels, rotation=30)

        fig.set_size_inches(15,10)
        fig.set_dpi(300)

#	fig.savefig(''.join([OUTDIR,'MEAN_Hov_201501_',st_wmo,'.png']))
	fig.savefig(''.join([OUTDIR,title,'.png']))
np.savetxt(''.join([OUTDIR,'AUG_persistence_index.txt']),table_integ,fmt="%10.5f",delimiter="\t",newline="\n")
np.savetxt(''.join([OUTDIR,'AUG_persistence_index2.txt']),table_integ2,fmt="%10.5f",delimiter="\t",newline="\n")
sys.exit()
