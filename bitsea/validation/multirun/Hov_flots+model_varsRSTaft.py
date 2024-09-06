import numpy as np
import numpy.ma as ma
import pylab as plt
import matplotlib.dates as mdates

from basins import OGS
from instruments.matchup_manager import Matchup_Manager
from timeseries.plot import *
from profilerRSTaft_2017 import *
from instruments import check
from instruments.var_conversions import FLOATVARS
# from instruments.var_conversions import LOVFLOATVARS as FLOATVARS
from instruments import superfloat as bio_float
# from instruments import lovbio_float as bio_float
from layer_integral import coastline
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
import scipy.io.netcdf as NC
import argparse
from commons.utils import addsep


def argument():
    parser = argparse.ArgumentParser(description='''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--maskfile', '-m',
                        type=str,
                        default="/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
                        required=False,
                        help=''' Path of maskfile''')

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help="")

    parser.add_argument('--statsdir', '-i',
                        type=str,
                        default=None,
                        required=False,
                        help="")

    return parser.parse_args()


args = argument()


TheMask = Mask(args.maskfile)
OUTDIR = addsep(args.outdir)

plotlines = False
if not (args.statsdir is None):
    DIRSTATS = addsep(args.statsdir)
    plotlines = True

Check_obj = check.check("", verboselevel=0)

font_s = 15
label_s = 15


def readModelProfile(filename, var, wmo):
    ncIN = NC.netcdf_file(filename, 'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile, :]
    return Profile


def get_level300(TheMask):

    i0 = TheMask.getDepthIndex(300)
    i1 = i0 + 1

    diff0 = 300 - TheMask.zlevels[i0]
    diff1 = 300 - TheMask.zlevels[i1]
    data = [(i0, diff0), (i1, diff1)]
    ix, datamin = min(data, key=lambda t: t[1])

    return ix


T_start = DATESTART
T_end = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start, '%Y%m%d'))
T_end2num = mpldates.date2num(datetime.strptime(T_end, '%Y%m%d'))
reg1 = [OGS.med]
reg_sn = ['med']

#max_depth = 26
# max_depth = 44 # 303m in the 142 levels model
# max_depth = get_level300(TheMask)

MM = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
#varname = ['CHLA','DOXY','NITRATE','TEMP','PSAL']
plotvarname = [r'Chl $[mg/m^3]$', r'Oxy $[mmol/m^3]$',
               r'Nitr $[mmol/m^3]$']  # ,r'Temp $[^\circ C]$','Sal']
#read_adjusted = [True,False,False,False,False]
#mapgraph = [3,4,5,1,2]
VARLIST = ['P_l', 'O2o', 'N3n']
VARLIST = ['Chla', 'O2o', 'N3n']
VARLIST = ['Chla', 'O2o', 'N3n']
Adj = [True, False, True]
VARLIST = ['P_l']
Adj = [True]
VARLIST = ['P_l', 'N3n']
Adj = [True, True]
VARLIST = ['N3n','P_l']
# Adj = [True]
nVar = len(VARLIST)


depthDICT = {
    'P_l': 300,
    'N3n': 500,
}
for ivar, var_mod in enumerate(VARLIST):
    var = FLOATVARS[var_mod]
    Profilelist_1 = bio_float.FloatSelector(var, TI1, OGS.med)
    wmo_list = bio_float.get_wmo_list(Profilelist_1)
    print '------------' + var_mod + '-----------'
    print wmo_list
    ndepths = depthDICT[var_mod]/5+1
    NewPres_5m = np.linspace(0, depthDICT[var_mod], ndepths)
    depths = NewPres_5m

    for wmo in wmo_list:
        if plotlines:
            filestats = DIRSTATS + '/' + wmo + '.nc'
            stats = ncreader(filestats)
        list_float_track = bio_float.filter_by_wmo(Profilelist_1, wmo)
        nP = len(list_float_track)
        Lon = np.zeros((nP,), np.float64)
        Lat = np.zeros((nP,), np.float64)

        plotmat = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_model = np.zeros([len(depths), len(list_float_track)])*np.nan
        timelabel_list = list()
        # adj = Adj[ivar]

        for ip, p in enumerate(list_float_track):
            Lon[ip] = p.lon
            Lat[ip] = p.lat
            Pres, Prof, Qc = p.read(var) #, read_adjusted=adj)
            # Pres, Prof, Qc = p.read(var, read_adjusted=adj)
            ii = Pres <= 300
            # Deve avere almeno 5 records:
            timelabel_list.append(p.time)
            if len(Prof[ii]) < 5:
                continue
                # NewProf_5m = np.interp(NewPres_5m, Pres[ii], Prof[ii])
                # plotmat[:, ip] = NewProf_5m

            try:
                GM = MM.getMatchups2([p], TheMask.zlevels, var_mod, interpolation_on_Float=False, \
                    checkobj=Check_obj, forced_depth=depths, extrapolation=False)
                # TM = MM.modeltime(p)
                # FILENAME = BASEDIR + \
                # TM.strftime("PROFILES/ave.%Y%m%d-12:00:00.profiles.nc")
                # M = readModelProfile(FILENAME, var_mod, p.ID())
                # M_newDepth = np.interp(
                #     NewPres_5m, TheMask.zlevels[:max_depth+1], M[:max_depth+1])
                plotmat_model[:, ip] = GM.Model
                plotmat[      :, ip] = GM.Ref
            except:
                print ' ... Not exists ' + p.time.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
                continue

        print var_mod + " " + np.str(len(timelabel_list)) + p.available_params
        if p.available_params.find(var) < 0:
            continue

        fig = plt.figure()

        ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2)
        ax2 = plt.subplot2grid((3, 3), (0, 2))
        ax3 = plt.subplot2grid((3, 3), (1, 0), colspan=3)
        ax4 = plt.subplot2grid((3, 3), (2, 0), colspan=3)

        xs, ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_m = ma.masked_invalid(plotmat)
        if (var_mod == 'P_l') or (var_mod == 'Chla'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat_m, shading='flat',
                                      vmin=0.00, vmax=0.40, cmap="viridis")  # default is 'flat'
            if plotlines:
                _, ref_dcm = stats.plotdata(var_mod, 'DCM')
                ax3.plot(mpldates.date2num(timelabel_list), ref_dcm, '-k')
                _, ref_mld = stats.plotdata(var_mod, 'MLD3')
                ax3.plot(mpldates.date2num(timelabel_list), ref_mld, '-w')
        if (var_mod == 'O2o'):
            quadmesh = ax3.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=200, vmax=250)
        if (var_mod == 'N3n'):
            quadmesh = ax3.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=0.00, vmax=6)
            if plotlines:
                _, ref_nit3 = stats.plotdata(var_mod, 'Nit_prof')
                ax3.plot(mpldates.date2num(timelabel_list), ref_nit3, '-k')
        # Inform matplotlib that the x axis is made by dates
        ax3.set_xlim([T_start2num, T_end2num])
        ax3.set_xticklabels([])
        ax3.set_ylim(0, depths[-1])
        ax3.invert_yaxis()
        ax3.set_ylabel("depth $[m]$", color='k', fontsize=font_s)
        fig.set_size_inches(15, 10)

        ###ax = axs[2]
        xs, ys = np.meshgrid(mpldates.date2num(timelabel_list), depths)
        plotmat_model_m = ma.masked_invalid(plotmat_model)
        if (var_mod == 'P_l') or (var_mod == 'Chla'):
            quadmesh = ax4.pcolormesh(xs, ys, plotmat_model_m, shading='flat',
                                      vmin=0.00, vmax=0.40, cmap="viridis")  # default is 'flat'
            if plotlines:
                model_dcm, _ = stats.plotdata(var_mod, 'DCM')
                ax4.plot(mpldates.date2num(timelabel_list), model_dcm, '-k')
                model_mld, _ = stats.plotdata(var_mod, 'MLD3')
                ax4.plot(mpldates.date2num(timelabel_list), model_mld, '-w')
        if (var_mod == 'O2o'):
            quadmesh = ax4.pcolormesh(
                xs, ys, plotmat_model_m, shading='flat', vmin=200, vmax=250)
        if (var_mod == 'N3n'):
            quadmesh = ax4.pcolormesh(
                xs, ys, plotmat_model_m, shading='flat', vmin=0.00, vmax=6)
            if plotlines:
                model_nit3, _ = stats.plotdata(var_mod, 'Nit_prof')
                ax4.plot(mpldates.date2num(timelabel_list), model_nit3, '-k')
        ax4.set_xlim([T_start2num, T_end2num])
        ax4.xaxis_date()
        ax4.set_ylim(0, depths[-1])
        ax4.invert_yaxis()
        ax4.set_ylabel("depth $[m]$", color='k', fontsize=font_s)
        ax4.xaxis.set_major_locator(
            mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11]))
        ax4.xaxis.set_major_formatter(
            mdates.DateFormatter("%Y%m"))  # ("%b%y"))
        for tick in ax4.get_xticklabels():
            tick.set_rotation(45)

        ax1.tick_params(axis='both', which='major', labelsize=label_s)
        ax2.tick_params(axis='both', which='major', labelsize=label_s)
        ax3.tick_params(axis='both', which='major', labelsize=label_s)
        ax4.tick_params(axis='both', which='major', labelsize=label_s)

        colorbar_bottom = ax4.get_position().ymin
        colorbar_extent = ax3.get_position().ymax - colorbar_bottom

        cbaxes = fig.add_axes([0.93, colorbar_bottom, 0.02, colorbar_extent])
        cbar = fig.colorbar(quadmesh, cax=cbaxes, )
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)

        c_lon, c_lat = coastline.get()

        ax1.plot(c_lon, c_lat, 'k')
        ax1.plot(Lon, Lat, 'r.')
        ax1.set_title("TRAJECTORY of FLOAT " + p.name() +
                      " - " + var_mod, color='r', fontsize=font_s)
        ind_max_sup = plotmat[0, :].argmax()
        ax1.plot(Lon[0], Lat[0], 'bo')
        ax1.plot(Lon[0], Lat[0], 'bx')
        ax1.set_xlim([-5.5, 36])
        extent = 4  # degrees
        ipp = len(list_float_track)
        ax2.plot(c_lon, c_lat, 'k')
        ax2.plot(Lon, Lat, 'ro')
        ax2.plot(Lon[0], Lat[0], 'bo')
        ax2.plot(Lon[0], Lat[0], 'bx', markersize=font_s)
        ax2.set_xlim([np.min(Lon[:ipp]) - extent/2,
                      np.max(Lon[:ipp]) + extent/2])
        ax2.set_ylim([np.min(Lat[:ipp]) - extent/2,
                      np.max(Lat[:ipp]) + extent/2])

        fig.savefig(
            ''.join([OUTDIR, 'HovRSTaft_Float+TRANS_', p.name(), '_', var_mod, '.png']))
        plt.close(fig)
