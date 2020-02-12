# , add_metadata
import sys
import os
import argparse
import numpy as np
import numpy.ma as ma
import instruments
import matplotlib.dates as mdates
import matplotlib.pyplot as pl
import scipy.io.netcdf as NC

from profiler_2015 import *

from commons.utils import addsep
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask

from instruments.matchup_manager import Matchup_Manager
from instruments.var_conversions import LOVFLOATVARS
from instruments import lovbio_float as bio_float
from instruments import matchup_manager
from layer_integral import coastline
from mhelpers.pgmean import PLGaussianMean
from timeseries.plot import *

from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader
from validation.online.profileplotter import figure_generator, ncwriter
from basins import OGS
# from basins.region import Region, Rectangle
# import matplotlib


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


meanObj11 = PLGaussianMean(5, 1.0)

TheMask = Mask(args.maskfile)
OUTDIR = addsep(args.outdir)

plotlines = False
if not (args.statsdir is None):
    plotlines = True
    DIRSTATS = addsep(args.statsdir)

font_s = 15
label_s = 15


def my_Hovmoeller_diagram(plotmat, xs, ys, fig=None, ax=None):
    if (fig is None) or (ax is None):
        fig, ax = pl.subplots()
    quadmesh = ax.pcolormesh(
        xs, ys, plotmat, shading='gouraud')  # default is 'flat'
    # Inform matplotlib that the x axis is made by dates
    ax.xaxis_date()
    ax.invert_yaxis()
    return fig, ax, quadmesh


def readModelProfile(filename, var, wmo):
    ncIN = NC.netcdf_file(filename, 'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile, :]
    return Profile


def get_leveldep(TheMask, dep):

    i0 = TheMask.getDepthIndex(dep)
    i1 = i0 + 1

    diff0 = dep - TheMask.zlevels[i0]
    diff1 = dep - TheMask.zlevels[i1]
    data = [(i0, diff0), (i1, diff1)]
    ix, datamin = min(data, key=lambda t: t[1])

    return ix


T_start = DATESTART
T_end = DATE__END
TI1 = T_INT
T_start2num = mpldates.date2num(datetime.strptime(T_start, '%Y%m%d'))
T_end2num = mpldates.date2num(datetime.strptime(T_end, '%Y%m%d'))
# reg1 = [OGS.med]
# reg_sn = ['med']

# max_depth = 26
# max_depth = 44 # 303m in the 142 levels model

MM = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
# varname = ['CHLA','DOXY','NITRATE','TEMP','PSAL']
plotvarname = [r'Chl $[mg/m^3]$', r'Oxy $[mmol/m^3]$',
               r'Nitr $[mmol/m^3]$']  # ,r'Temp $[^\circ C]$','Sal']
# mapgraph = [3,4,5,1,2]
VARLIST = ['P_l', 'O2o', 'N3n']
VARLIST = ['Chla', 'O2o', 'N3n']
VARLIST = ['N3n']

DICTadj = {
    'P_l': True,
    'Chla': True,
    'O2o': False,
    'N3n': True,
}


nVar = len(VARLIST)

meanObj11 = PLGaussianMean(11, 1.0)

Profilelist_1 = bio_float.FloatSelector(None, TI1, OGS.med)
wmo_list = bio_float.get_wmo_list(Profilelist_1)

print wmo_list

# NewPres_5m=np.linspace(0,300,61)
NewPres_5m = np.linspace(0, 40, 9)
NewPres_10m = np.linspace(50, 200, (200-50)/10+1)
NewPres_25m = np.linspace(225, 1000, (1000-225)/25+1)
NewPres = np.concatenate([NewPres_5m, NewPres_10m, NewPres_25m])

dep1 = 300
dep2 = 650
maskPres_d1 = NewPres <= dep1
maskPres_d2 = (NewPres >= dep1) & (NewPres <= dep2)
depths1 = NewPres[maskPres_d1]
depths2 = NewPres[maskPres_d2]

dep_plot = 650
maskPresplot = NewPres <= dep_plot

depths = NewPres[maskPresplot]
maskDepth1 = depths <= dep1
maskDepth2 = (depths >= dep1) & (depths<= dep2)

LISTlayers = [
    '200m',
    '200_400m',
    '400_600m',
]

maskPres = {
    '200m': NewPres <= 200,
    '200_400m': (NewPres > 200) & (NewPres <= 400),
    '400_600m': (NewPres > 400) & (NewPres <= 600),
}

fixedThreshold = {
    '200m': 3,
    '200_400m': 4,
    '400_600m': 5,
}

#nL = len(LISTlayers)



# for j in range(0,len(wmo_list)):
LIST = [[] for ii in range(5)]
acceptN3n = {}
for j, wmo in enumerate(wmo_list):
    # if not(wmo[-2]=='7'): continue
    print wmo
    #    if (j==0):
    if plotlines:
        filestats = DIRSTATS + '/' + wmo + '.nc'
        stats = ncreader(filestats)
    list_float_track = bio_float.filter_by_wmo(Profilelist_1, wmo)
    nP = len(list_float_track)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)

    acceptN3n[wmo] = np.zeros(nP)
    acceptN3n[wmo][:] = np.nan
    for ivar, var_mod in enumerate(VARLIST):
        print ' ... ' + var_mod
        #       if (ivar==0):
        plotmat = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_model = np.zeros([len(depths), len(list_float_track)])*np.nan
        plotmat_misf = np.zeros([len(depths), len(list_float_track)])*np.nan
        timelabel_list = list()
        var = LOVFLOATVARS[var_mod]
        adj = DICTadj[var_mod]

        for ip, p in enumerate(list_float_track):
            Lon[ip] = p.lon
            Lat[ip] = p.lat
            Pres, Prof, Qc = p.read(var, read_adjusted=adj)
            ii = Pres <= 300
# Deve avere almeno 5 records:
            if len(Prof[ii]) > 5:
                NewProf = np.interp(NewPres, Pres, Prof)
                plotmat[:, ip] = NewProf[maskPresplot]
# FILTRAGGIO per CHLA e N3n dato dai valori troppo alti registrato in superficie:
                if (var_mod == "P_l"):
                    if (plotmat[0, ip] > 0.45):
                        plotmat[:, ip] = np.nan
                if (var_mod == "N3n"):
                    if (plotmat[0, ip] > 2):
                        plotmat[:, ip] = np.nan
                        NewProf[:] = np.nan
                    else:
                        acceptN3n[wmo][ip] = 1


            timelabel_list.append(p.time)

            # PLOT FOR THE MODEL
            TM = MM.modeltime(p)
            FILENAME = BASEDIR + \
                TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
            M = readModelProfile(FILENAME, var_mod, p.ID())
            M_newDepth = np.interp(NewPres, TheMask.zlevels, M)
            plotmat_model[:, ip] = M_newDepth[maskPresplot]

            plotmat_misf[:,ip] = plotmat[:,ip]-plotmat_model[:,ip]


            # Blacklisting
            if (var_mod == 'N3n') and (acceptN3n[wmo][ip] == 1):
                misf2 = (NewProf - M_newDepth)**2
                mean200 = np.nanmean(NewProf[maskPres['200m']])
                for ll,lay in enumerate(LISTlayers):
                    maskPres[lay]
                    rmsd = (np.nanmean(misf2[maskPres[lay]]))**2
                    if rmsd/mean200 > fixedThreshold[lay]:
                        plotmat[:,ip] = np.nan
                  
        if np.all(np.isnan(plotmat_misf)): continue

        print var_mod + " " + np.str(len(timelabel_list)) + p.available_params

        fig = pl.figure()
        fig.set_size_inches(15,10)

        ax1 = pl.subplot2grid((3, 3), (0, 0), colspan=2)
        ax2 = pl.subplot2grid((3, 3), (0, 2))
        ax3 = pl.subplot2grid((3, 3), (1, 0), colspan=3)
        ax4 = pl.subplot2grid((3, 3), (2, 0), colspan=3)

        xs, ys = np.meshgrid(mpldates.date2num(timelabel_list), depths1)
        plotmat_m = ma.masked_invalid(plotmat[maskDepth1,:])
        if (var_mod == 'P_l') or (var_mod == 'Chla'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat_m, shading='flat',
                                      vmin=0.00, vmax=0.40, cmap="viridis")  # default is 'flat'
            if plotlines:
                _, ref_dcm = stats.plotdata(var_mod, 'DCM')
                ax3.plot(mpldates.date2num(timelabel_list), ref_dcm, '-k')
                _, ref_mld = stats.plotdata(var_mod, 'MLD3')
                ax3.plot(mpldates.date2num(timelabel_list), ref_mld, '-w')
        if (var_mod == 'O2o'):
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax3.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=200, vmax=250)
        if (var_mod == 'N3n'):
            # print len(xs)
            # print len(ys)
            # print plotmat_m.shape
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax3.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=0.00, vmax=4)
            if plotlines:
                _, ref_nit3 = stats.plotdata(var_mod, 'Nit_prof')
                ax3.plot(mpldates.date2num(timelabel_list), ref_nit3, '-k')
        # Inform matplotlib that the x axis is made by dates
        ax3.set_xlim([T_start2num, T_end2num])
#	ax3.xaxis_date()
        ax3.set_xticklabels([])
        ax3.set_ylim(0, dep1)
        ax3.invert_yaxis()
#	ax3.set_title("FLOAT " + p.name() + ": " + plotvarname[ivar], color = 'b')
        ax3.set_ylabel("depth $[m]$", color='k', fontsize=font_s)
#	fig.set_dpi(300)

        # ax = axs[2]
        xs, ys = np.meshgrid(mpldates.date2num(timelabel_list), depths2)
        plotmat_m = ma.masked_invalid(plotmat[maskDepth2,:])
        if (var_mod == 'P_l') or (var_mod == 'Chla'):
            quadmesh = ax4.pcolormesh(xs, ys, plotmat_m, shading='flat',
                                      vmin=0.00, vmax=0.40, cmap="viridis")  # default is 'flat'
            if plotlines:
                _, ref_dcm = stats.plotdata(var_mod, 'DCM')
                ax3.plot(mpldates.date2num(timelabel_list), ref_dcm, '-k')
                _, ref_mld = stats.plotdata(var_mod, 'MLD3')
                ax3.plot(mpldates.date2num(timelabel_list), ref_mld, '-w')
        if (var_mod == 'O2o'):
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax4.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=200, vmax=250)
        if (var_mod == 'N3n'):
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax4.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=0.00, vmax=4)
            if plotlines:
                _, ref_nit3 = stats.plotdata(var_mod, 'Nit_prof')
                ax3.plot(mpldates.date2num(timelabel_list), ref_nit3, '-k')
        ax4.set_xlim([T_start2num, T_end2num])
        ax4.xaxis_date()
        ax4.set_ylim(dep1, dep2)
        ax4.invert_yaxis()
#	ax4.set_title("BFM model: " + plotvarname[ivar], color = 'b')
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

  #      cbaxes = fig.add_axes([0.95, 0.01, 0.02, 0.95])
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
#	ax1.set_ylabel("LAT",color = 'k',fontsize=font_s)
#	ax1.set_xlabel("LON",color = 'k')
        extent = 4  # degrees
        ipp = len(list_float_track)
        ax2.plot(c_lon, c_lat, 'k')
        ax2.plot(Lon, Lat, 'ro')
#	ax2.plot(Lon[ind_max_sup],Lat[ind_max_sup],'go')
        ax2.plot(Lon[0], Lat[0], 'bo')
        ax2.plot(Lon[0], Lat[0], 'bx', markersize=font_s)
        ax2.set_xlim([np.min(Lon[:ipp]) - extent/2,
                      np.max(Lon[:ipp]) + extent/2])
        ax2.set_ylim([np.min(Lat[:ipp]) - extent/2,
                      np.max(Lat[:ipp]) + extent/2])

        fig.savefig(
            ''.join([OUTDIR, 'Hov_Floatexcl_', p.name(), '_', var_mod, '.png']))
        pl.close(fig)



        fig = pl.figure()
        fig.set_size_inches(15,10)

        ax1 = pl.subplot2grid((3, 3), (0, 0), colspan=2)
        ax2 = pl.subplot2grid((3, 3), (0, 2))
        ax3 = pl.subplot2grid((3, 3), (1, 0), colspan=3)
        ax4 = pl.subplot2grid((3, 3), (2, 0), colspan=3)

        xs, ys = np.meshgrid(mpldates.date2num(timelabel_list), depths1)
        plotmat_m = ma.masked_invalid(plotmat_misf[maskDepth1,:])
        if (var_mod == 'P_l') or (var_mod == 'Chla'):
            quadmesh = ax3.pcolormesh(xs, ys, plotmat_m, shading='flat',
                                      vmin=0.00, vmax=0.40, cmap="viridis")  # default is 'flat'
            if plotlines:
                _, ref_dcm = stats.plotdata(var_mod, 'DCM')
                ax3.plot(mpldates.date2num(timelabel_list), ref_dcm, '-k')
                _, ref_mld = stats.plotdata(var_mod, 'MLD3')
                ax3.plot(mpldates.date2num(timelabel_list), ref_mld, '-w')
        if (var_mod == 'O2o'):
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax3.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=200, vmax=250)
        if (var_mod == 'N3n'):
            # print len(xs)
            # print len(ys)
            # print plotmat_m.shape
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax3.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=0.00, vmax=4)
            if plotlines:
                _, ref_nit3 = stats.plotdata(var_mod, 'Nit_prof')
                ax3.plot(mpldates.date2num(timelabel_list), ref_nit3, '-k')
        # Inform matplotlib that the x axis is made by dates
        ax3.set_xlim([T_start2num, T_end2num])
#	ax3.xaxis_date()
        ax3.set_xticklabels([])
        ax3.set_ylim(0, dep1)
        ax3.invert_yaxis()
#	ax3.set_title("FLOAT " + p.name() + ": " + plotvarname[ivar], color = 'b')
        ax3.set_ylabel("depth $[m]$", color='k', fontsize=font_s)
#	fig.set_dpi(300)

        # ax = axs[2]
        xs, ys = np.meshgrid(mpldates.date2num(timelabel_list), depths2)
        plotmat_m = ma.masked_invalid(plotmat_misf[maskDepth2,:])
        if (var_mod == 'P_l') or (var_mod == 'Chla'):
            quadmesh = ax4.pcolormesh(xs, ys, plotmat_m, shading='flat',
                                      vmin=0.00, vmax=0.40, cmap="viridis")  # default is 'flat'
            if plotlines:
                _, ref_dcm = stats.plotdata(var_mod, 'DCM')
                ax3.plot(mpldates.date2num(timelabel_list), ref_dcm, '-k')
                _, ref_mld = stats.plotdata(var_mod, 'MLD3')
                ax3.plot(mpldates.date2num(timelabel_list), ref_mld, '-w')
        if (var_mod == 'O2o'):
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax4.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=200, vmax=250)
        if (var_mod == 'N3n'):
            # ,cmap="jet")# default is 'flat'
            quadmesh = ax4.pcolormesh(
                xs, ys, plotmat_m, shading='flat', vmin=0.00, vmax=4)
            if plotlines:
                _, ref_nit3 = stats.plotdata(var_mod, 'Nit_prof')
                ax3.plot(mpldates.date2num(timelabel_list), ref_nit3, '-k')
        ax4.set_xlim([T_start2num, T_end2num])
        ax4.xaxis_date()
        ax4.set_ylim(dep1, dep2)
        ax4.invert_yaxis()
#	ax4.set_title("BFM model: " + plotvarname[ivar], color = 'b')
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

  #      cbaxes = fig.add_axes([0.95, 0.01, 0.02, 0.95])
        cbaxes = fig.add_axes([0.93, colorbar_bottom, 0.02, colorbar_extent])
        cbar = fig.colorbar(quadmesh, cax=cbaxes, )
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=font_s)

        c_lon, c_lat = coastline.get()

        ax1.plot(c_lon, c_lat, 'k')
        ax1.plot(Lon, Lat, 'r.')
        ax1.set_title("TRAJECTORY of FLOAT " + p.name() +
                      " - " + var_mod, color='r', fontsize=font_s)
        ind_max_sup = plotmat_misf[0, :].argmax()
        ax1.plot(Lon[0], Lat[0], 'bo')
        ax1.plot(Lon[0], Lat[0], 'bx')
        ax1.set_xlim([-5.5, 36])
#	ax1.set_ylabel("LAT",color = 'k',fontsize=font_s)
#	ax1.set_xlabel("LON",color = 'k')
        extent = 4  # degrees
        ipp = len(list_float_track)
        ax2.plot(c_lon, c_lat, 'k')
        ax2.plot(Lon, Lat, 'ro')
#	ax2.plot(Lon[ind_max_sup],Lat[ind_max_sup],'go')
        ax2.plot(Lon[0], Lat[0], 'bo')
        ax2.plot(Lon[0], Lat[0], 'bx', markersize=font_s)
        ax2.set_xlim([np.min(Lon[:ipp]) - extent/2,
                      np.max(Lon[:ipp]) + extent/2])
        ax2.set_ylim([np.min(Lat[:ipp]) - extent/2,
                      np.max(Lat[:ipp]) + extent/2])

        fig.savefig(
            ''.join([OUTDIR, 'Hov_misf_', p.name(), '_', var_mod, '.png']))
        pl.close(fig)



# Blacklisting plots


# LISTwmo = np.unique(LIST[0])



# meanNorm = np.zeros((nL, len(LISTwmo)))
# meanNorm[:,:] = np.nan


# thresholds = {
#     0: [aa/10. for aa in range(20,4,-5)],
#     1: [aa/10. for aa in range(30,14,-5)],
#     2: [aa/10. for aa in range(35,19,-5)],
# }

# nexc = {}
# Np = {}
# for wmo in LISTwmo:
#     maskwmo = np.array(LIST[0]) == wmo
#     Np[wmo] = np.sum(maskwmo)
#     nexc[wmo] = {}
#     for ii in range(nL):
#         nexc[wmo][ii] = np.zeros(len(thresholds[ii]))
#         nexc[wmo][ii][:] = np.nan
#         values = np.array(LIST[ii+1])[maskwmo]
#         for iit,tt in enumerate(thresholds[ii][1:]):
#             H_tt = thresholds[ii][iit]
#             nexc[wmo][ii][iit+1] = np.nansum((values>tt) & (values<H_tt))
#         nexc[wmo][ii][0] = np.nansum(values>thresholds[ii][0])



# thresholdsN = {
#     0: [aa/10. for aa in range(40,9,-5)],
#     1: [aa/10. for aa in range(50,19,-5)],
#     2: [aa/10. for aa in range(50,19,-5)],
# }


# nexcN = {}
# excprof = {}
# for wmo in LISTwmo:
#     maskwmo = np.array(LIST[0]) == wmo
#     nexcN[wmo] = {}
#     excprof[wmo] = np.zeros(Np[wmo])
#     for ii in range(nL):
#         nexcN[wmo][ii] = np.zeros(len(thresholdsN[ii]))
#         nexcN[wmo][ii][:] = np.nan
#         values = (np.array(LIST[ii+1])/np.array(LIST[4]))[maskwmo]
#         for iit,tt in enumerate(thresholdsN[ii][1:]):
#             H_tt = thresholdsN[ii][iit]
#             nexcN[wmo][ii][iit+1] = np.nansum((values>tt) & (values<H_tt))
#         nexcN[wmo][ii][0] = np.nansum(values>thresholdsN[ii][0])

#         mask_fixt = values>fixedThreshold[ii]
#         excprof[wmo][mask_fixt] += 2**ii 



# percexc = np.zeros(len(LISTwmo))
# percexc[:] = np.nan
# for iw,wmo in enumerate(LISTwmo):
#     percexc[iw] = 100.*np.nansum(excprof[wmo]>0)/Np[wmo]


# wlayerDICT= {
#     1: '200',
#     2: '2-400',
#     4: '4-600',
#     3: '200 2-400',
#     6: '2-400 4-600',
#     5: '200 4-600',
#     7: 'all',
# }

# layersH = np.zeros((len(LISTwmo),7))
# layersH[:,:] = np.nan
# layersB = np.zeros((len(LISTwmo),7))
# layersB[:,:] = np.nan
# layersB[:,0] = 0.
# for iw,wmo in enumerate(LISTwmo):
#     for ir in wlayerDICT.keys():
#         if np.nansum(excprof[wmo]>0)==0: continue
#         layersH[iw,ir-1] = 100.*np.nansum(excprof[wmo]==ir)/np.nansum(excprof[wmo]>0)
#     for irr in range(1,len(wlayerDICT.keys())):
#         layersB[iw,irr] = np.nansum(layersH[iw,:irr])




