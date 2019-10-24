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
    parser.add_argument(   '--basedir', '-b',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
                                required = True,
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



#import matplotlib
#matplotlib.use('Agg')
import os,sys
import scipy.io.netcdf as NC
from commons.mask import Mask
from commons.Timelist import TimeList, TimeInterval
from basins.region import Rectangle
from layer_integral import coastline
from instruments import superfloat as bio_float
from instruments.var_conversions import FLOATVARS
import numpy as np
import numpy.ma as ma
from commons.utils import addsep
import pylab as pl
from mhelpers.pgmean import PLGaussianMean
from instruments.matchup_manager import Matchup_Manager
import matplotlib.dates as mdates
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
import datetime
from mpl_toolkits.basemap import Basemap

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


xlim=[-6.5,36.5]
ylim=[30,46]
xC=(xlim[0]+xlim[1])/2
yC=(ylim[0]+ylim[1])/2
 
mapobj = Basemap(projection='merc',lat_0=xC,lon_0=yC,\
                                  llcrnrlon = xlim[0], \
                                  llcrnrlat = ylim[0], \
                                  urcrnrlon = xlim[1], \
                                  urcrnrlat = ylim[1], \
                                  area_thresh=None, \
                                  resolution='l')
background_color=(.9, .9, .9)
POLY=mapobj.coastpolygons
R1 = Rectangle(28.90338, 34.88078, 37.28154, 39.795254)#anatolia
R2 = Rectangle(7.86973,11.0276,45.1554,46 )#nord italia
R3 = Rectangle(11.7043,12.3246,42.3872,43.5016)# centro italia
R4 = Rectangle(-0.532455,0.369799,40.9977,41.8014)# catalunya
R5 = Rectangle(-6,-3.57754,39.011,42.0531)# centro spagna
R6 = Rectangle(20.332,21.6854,40.3991,41.4642) #albania
R7 = Rectangle(35.1627,35.9522,30.7444,33.0422)# israel
R8 = Rectangle(28.2831,28.9034,45.1951,45.7487)# moldavia
R9 = Rectangle(19.0258,19.588,42.0524,42.3803)# Skadarko Jezero
R10 = Rectangle(27.7687,28.845,40.0348,40.2938)# 3 laghi vicino bursa, turchia
R11 = Rectangle(29.2417,30.0876,39.9683,40.6135)# terzo lago
R12 = Rectangle(31.7135,32.3184,31.0228,31.5918)#nilo
RECTANGLE_LIST=[R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11] # R12]
def is_in_boxes(x,y,Rectangle_list, map_obj):
    for rect in Rectangle_list:
        lonmin,latmin =map_obj(rect.lonmin, rect.latmin)
        lonmax,latmax =map_obj(rect.lonmax, rect.latmax)
        rect_meters = Rectangle(lonmin, lonmax, latmin, latmax)
        if rect_meters.is_inside(x,y):
            return True
    return False


INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)
BASEDIR = addsep(args.basedir)
TheMask=Mask(args.maskfile)

font_s =  13
font_s2 = 3
label_s = 15


#######
# For second plot:




TL = TimeList.fromfilenames(None, BASEDIR + "PROFILES/","ave*.nc")
deltaT= datetime.timedelta(hours=12)
TI = TimeInterval.fromdatetimes(TL.Timelist[0] - deltaT, TL.Timelist[-1] + deltaT)
ALL_PROFILES = bio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))

METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1','SurfVal','dNit_dz']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
wmo_list=bio_float.get_wmo_list(ALL_PROFILES)
izmax = TheMask.getDepthIndex(200) 

################
T_start2num = mdates.date2num(TI.start_time)
T_end2num   = mdates.date2num(TI.end_time)

max_depth = get_level_depth(TheMask,300)

plotvarname = [r'Chl $[mg/m^3]$',r'Oxy $[mmol/m^3]$',r'Nitr $[mmol/m^3]$'] 
VARLIST = ['P_l','O2o','N3n']

nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11,1.0)

def var_counter(profilelist, var):
    counter=0
    for p in profilelist:
        if var in p.available_params:
            counter+=1
    return counter

bt=300
NewPres_5m=np.linspace(0,300,121)
depths=NewPres_5m
for ivar, var_mod in enumerate(VARLIST):
    var = FLOATVARS[var_mod]
    Profilelist = bio_float.FloatSelector(var, TI, Rectangle(-6,36,30,46))
    wmo_list=bio_float.get_wmo_list(Profilelist)
    
    for wmo in wmo_list:
        list_float_track=bio_float.filter_by_wmo(Profilelist,wmo)
        nP = len(list_float_track)
        if nP <2 : continue
        
        Lon = np.zeros((nP,), np.float64)
        Lat = np.zeros((nP,), np.float64)
        plotmat = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_model = np.zeros([len(depths), len(list_float_track)])*np.nan
        timelabel_list = list()
        
        
    	for ip, p in enumerate(list_float_track):
    	    Lon[ip] = p.lon
    	    Lat[ip] = p.lat
            try:
    	        Pres,Prof,Qc=p.read(var)
            except:
                timelabel_list.append(p.time)
                continue
    	    ii = Pres<=bt
            if ii.sum() < 2 :
                timelabel_list.append(p.time)
                continue
            NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])
            plotmat[:,ip]=NewProf_5m
            timelabel_list.append(p.time)

	    # PLOT FOR THE MODEL
            try:
                TM=M.modeltime(p)
                FILENAME = BASEDIR + TM.strftime("PROFILES/ave.%Y%m%d-12:00:00.profiles.nc")
                Modelprofile = M.readModelProfile(FILENAME,var_mod,p.ID())
            except:
                continue
            M_newDepth=np.interp(NewPres_5m,TheMask.zlevels[:max_depth+1],Modelprofile[:max_depth+1])
            plotmat_model[:,ip] = M_newDepth
            
        print var_mod + " " + np.str(len(timelabel_list)) +  p.available_params


        fig = pl.figure()
        fig.set_size_inches(12,15)
        fig.set_dpi(300)       
            # AXES FOR MAPS
        ax1 = pl.subplot2grid((7, 3), (0, 0), colspan=2)
        ax2 = pl.subplot2grid((7, 3), (0, 2))
            # AXES PLOT THE HOVMOELLER
        ax3 = pl.subplot2grid((7, 3), (1, 0), colspan=3)
        ax4 = pl.subplot2grid((7, 3), (2, 0), colspan=3)
        # AXES FOR STATISTICS
        ax5 = pl.subplot2grid((7, 3), (3, 0), colspan=3)
        ax6 = pl.subplot2grid((7, 3), (4, 0), colspan=3)
        ax7 = pl.subplot2grid((7, 3), (5, 0), colspan=3)
        ax8 = pl.subplot2grid((7, 3), (6, 0), colspan=3)

        # MAP OF TRAJECTORY:

#         for polygon in POLY:
#             x,y=polygon
#             ax1.fill(x,y,color=background_color)
#             if not is_in_boxes(x[0], y[0], RECTANGLE_LIST, mapobj):
#                 ax1.plot(x,y,linewidth=0.7, color="0.4")
        ax1.set_position([0.125, 0.786097560976,0.5, 0.16 ])   
        c_lon,c_lat=coastline.get()
        ax1.plot(c_lon,c_lat,'k')
        ax1.plot(Lon,Lat,'r.',markersize=font_s2)
        ax1.set_title("TRAJECTORY of FLOAT " + p.name() , color = 'r', fontsize=font_s)
        ind_max_sup=plotmat[0,:].argmax()
        ax1.plot(Lon[0],Lat[0],'bo',markersize=font_s2)
        ax1.plot(Lon[0],Lat[0],'bx')
        ax1.set_xlim([-5.5,36])
        extent=4 #degrees
        ax2.plot(c_lon,c_lat,'k')
        ax2.plot(Lon,Lat,'ro',markersize=font_s2)
        ax2.plot(Lon[0],Lat[0],'bo',markersize=font_s2)
        ax2.plot(Lon[0],Lat[0],'bx',markersize=font_s)
        ax2.set_xlim([np.min(Lon[:nP]) -extent/2, np.max(Lon[:nP]) +extent/2])
        ax2.set_ylim([np.min(Lat[:nP]) -extent/2, np.max(Lat[:nP]) +extent/2])



        xs,ys = np.meshgrid(mdates.date2num(timelabel_list), depths)
        plotmat_m=ma.masked_invalid(plotmat)     
        
        # PLOT HOVMOELLER OF FLOAT
        if (var_mod == 'P_l'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=0.00,vmax=0.40,cmap="viridis")# default is 'flat'
            quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=0.00,vmax=0.40,cmap="viridis")# default is 'flat'
        if (var_mod == 'O2o'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=160,vmax=250) #,cmap="jet")# default is 'flat'
            quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=160,vmax=250)
        if (var_mod == 'N3n'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat_m,shading='flat',vmin=0.00,vmax=4) #,cmap="jet")# default is 'flat'
            quadmesh = ax4.pcolormesh(xs, ys, plotmat_model,shading='flat',vmin=0.00,vmax=4) #,cmap="jet")# default is 'flat'
	#Inform matplotlib that the x axis is made by dates
        ax3.set_xlim([T_start2num,timelabel_list[-1]])  # SET THE END BY LENGHT OF TIMESERIES
        ax3.set_ylabel("FLOAT \n depth $[m]$",color = 'k', fontsize=font_s)        
        ax3.set_xlim([T_start2num,timelabel_list[-1]])
        ax3.invert_yaxis()
        ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))

        ax4.set_xlim([T_start2num,timelabel_list[-1]])
        ax4.invert_yaxis()
        ax4.set_ylabel("MODEL \ndepth $[m]$",color = 'k',fontsize=font_s)
        ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax4.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m")) #("%b%y"))
        #for tick in ax4.get_xticklabels(): tick.set_rotation(45)
        
        ax1.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax2.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax3.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax4.tick_params(axis = 'both', which = 'major', labelsize = label_s)

        colorbar_bottom= ax4.get_position().ymin
        colorbar_extent= ax3.get_position().ymax - colorbar_bottom

        cbaxes = fig.add_axes([0.93, colorbar_bottom , 0.02, colorbar_extent])

# SECOND FIGURE ON STATISTICS:

        # READ STATISTICS FILE *.nc
        INPUT_FILE = "%s%s_%s.nc" %(INDIR, var_mod, wmo )
        A = ncreader(INPUT_FILE)
        times = timelabel_list
        OUTFILE = OUTDIR + var_mod + "_" + wmo + ".png"
        print OUTFILE

        model, ref =           A.plotdata(var_mod,'Int_0-200', only_good=False)
        if np.isnan(ref).all(): continue
        surf_model, surf_ref = A.plotdata(var_mod,'SurfVal', only_good=False)
        model_corr , ref_corr =A.plotdata(var_mod,'Corr', only_good=False)



        ax6.plot(times,  ref,'.b',label='REF INTEG')
        ax6.plot(times,model,'b',label='MOD INTEG')
	    
        ax5.plot(times,  surf_ref,'.b',label='FLOAT') #'REF SURF')
        ax5.plot(times,  surf_model,'b',label='MODEL') #'MOD SURF')
        legend = ax5.legend(loc='upper left', shadow=True, fontsize=12)
        if ( var_mod == "P_l" ):
            ax6.set_ylabel('INTG 0-200 \n $[mg{\  } m^{-3}]$',fontsize=15)
            ax5.set_ylabel('SURF \n $[mg{\  } m^{-3}]$',fontsize=15)
            ax5.set_ylim(0,0.5)
            ax6.set_ylim(0,0.22)
            model_dcm, ref_dcm =A.plotdata(var_mod,'DCM', only_good=False)
            model_mld, ref_mld =A.plotdata(var_mod,'z_01',only_good=False)
            ax8.invert_yaxis()
            ax8.plot(times,  ref_dcm,'.b',label='DCM REF')
            ax8.plot(times,model_dcm,'b',label='DCM MOD')
            ax8.plot(times, ref_mld,'.r',label='WLB REF') # WINTER AYER BLOOM
            ax8.plot(times,model_mld,'r',label='WLB MOD')
            ax8.set_ylabel('DCM/WLB $[m]$',fontsize=15)
            #ax8.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m-%Y"))
            ax8.set_xlim([T_start2num,timelabel_list[-1]])
            ax8.set_ylim([200,0])
         


        if ( var_mod == "O2o" ):
            ax6.set_ylabel('INTG 0-200m \n $[mmol{\  } m^{-3}]$',fontsize=15)
            ax5.set_ylabel('SURF \n $[mmol{\  } m^{-3}]$',fontsize=15)
            ax8.plot(times,  np.ones_like(times))

                        
        if ( var_mod == "N3n" ):
            ax6.set_ylabel('INTG 0-350m \n $[mmol{\  } m^{-3}]$',fontsize=15)
            ax5.set_ylim(0,1.8)
            ax6.set_ylim(0,2.8)
            ax5.set_ylabel('SURF VAL.\n $[mmol{\  } m^{-3}]$',fontsize=15)
            ax8.set_ylim(0,350) # NITRICLINE RANGE
            model_nit, ref_nit =A.plotdata(var_mod,'Nit_1',only_good=False)
            model_nit2, ref_nit2 = A.plotdata(var_mod,'dNit_dz', only_good=False)
            ref_nit2[ np.isnan(ref_nit2) ] = 0.0
            jj=[ref_nit2[ ~np.isnan(ref_nit2) ] <= 10.0 ] 
            ref_nit2[tuple(jj)] = np.nan
            if (~np.isnan(ref_nit).all() == True) or (~np.isnan(model_nit).all() == True):
                ax8.plot(times,  ref_nit,'.b',label='REF')
                ax8.plot(times,model_nit,'b',label='MOD')
                ax8.invert_yaxis()
                ax8.set_ylabel('NITRCL 1/2 $[m]$',fontsize=15)
                ax8.plot(times, ref_nit2,'.r',label='dNit REF')
                ax8.plot(times,model_nit2,'r',label='dNit MOD') 
            else: 
                ax8.plot(times,  np.ones_like(times))
                
               
        
        times_r = times
        try:
            ax2.set_xticklabels([])
            ax3.set_xticklabels([])
            ax4.set_xticklabels([])
            ax5.set_xticklabels([])
            ax6.set_xticklabels([])
        except:
            print "nans in figure"
        if (np.isnan(ref_corr).all() == False ):
            ax7.plot(times,ref_corr,'b')
            ax7.set_ylabel('CORR',fontsize=15)
            ax7.set_xticklabels([])
            ax7.set_ylim(0,1)
            
            ax7.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
            ax7.tick_params(axis = 'both', which = 'major', labelsize = label_s)
            ax7.set_xlim([T_start2num,timelabel_list[-1]])



        cbar=fig.colorbar(quadmesh, cax = cbaxes)
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)

        ax5.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax6.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax8.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
        ax8.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
        xlabels = ax8.get_xticklabels()
        pl.setp(xlabels, rotation=30, fontsize=15)
        
        ax3.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax3.set_xlim([T_start2num,timelabel_list[-1]])
        ax4.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax4.set_xlim([T_start2num,timelabel_list[-1]])

        ax5.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax5.set_xlim([T_start2num,timelabel_list[-1]])
        ax6.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax6.set_xlim([T_start2num,timelabel_list[-1]])
        ax8.tick_params(axis = 'both', which = 'major', labelsize = label_s)
        ax8.set_xlim([T_start2num,timelabel_list[-1]])
        ax5.xaxis.grid(True)
        ax6.xaxis.grid(True)
        ax7.xaxis.grid(True)
        ax8.xaxis.grid(True)

        fig.savefig(OUTFILE)
        pl.close(fig)
        import sys
        sys.exit()
