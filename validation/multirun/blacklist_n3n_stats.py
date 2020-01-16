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

# max_depth = get_leveldep(TheMask, 300)
# max_depth600 = get_leveldep(TheMask, 600)

MM = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
plotvarname = [r'Chl $[mg/m^3]$', r'Oxy $[mmol/m^3]$',
               r'Nitr $[mmol/m^3]$']  # ,r'Temp $[^\circ C]$','Sal']
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
maskPres200 = NewPres <= 200
maskPres2_400 = (NewPres > 200) & (NewPres <= 400)
maskPres4_600 = (NewPres > 400) & (NewPres <= 600)

dep_plot = 600
maskPresplot = NewPres <= dep_plot
depths = NewPres[maskPresplot]

# for j in range(0,len(wmo_list)):
LIST = [[] for ii in range(5)]
acceptN3n = {}
for j, wmo in enumerate(wmo_list):
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
# FILTRAGGIO per N3n dato dai valori troppo alti registrato in superficie:
                if (var_mod == "N3n"):
                    if (NewProf[0] <= 2):
                        acceptN3n[wmo][ip] = 1
                    else:
                        NewProf[:] = np.nan

            timelabel_list.append(p.time)

            # MODEL
            TM = MM.modeltime(p)
            FILENAME = BASEDIR + \
                TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
            M = readModelProfile(FILENAME, var_mod, p.ID())
            M_newDepth = np.interp(NewPres, TheMask.zlevels, M)

            if (var_mod == 'N3n') and (acceptN3n[wmo][ip] == 1):
                LIST[0].append(wmo[-3:])

                misf2 = (NewProf - M_newDepth)**2
                rmsd200 = (np.nanmean(misf2[maskPres200]))**.5
                LIST[1].append(rmsd200)
                rmsd2_400 = (np.nanmean(misf2[maskPres2_400]))**.5
                LIST[2].append(rmsd2_400)
                rmsd4_600 = (np.nanmean(misf2[maskPres4_600]))**.5
                LIST[3].append(rmsd4_600)

                mean200 = np.nanmean(NewProf[maskPres200])
                LIST[4].append(mean200)


        print var_mod + " " + np.str(len(timelabel_list)) + p.available_params
        
    if acceptN3n[wmo][ip]==1: break

import sys
sys.exit(0)


# Blacklisting statistic plots

pl.close('all')

fig1,axs = pl.subplots(nL, 1, sharex=True, sharey=True, figsize=[12, 8])

LISTwmo = np.unique(LIST[0])

titleDICT = {
    0: '200m',
    1: '200-400m',
    2: '400-600m',
}

nL = len(titleDICT.keys())

mean = np.zeros((nL, len(LISTwmo)))
std = np.zeros((nL, len(LISTwmo)))
mean[:,:] = np.nan
std[:,:] = np.nan
for iw, wmo in enumerate(LISTwmo):
    maskwmo = np.array(LIST[0]) == wmo
    for ii in range(nL):
        mean[ii, iw] = np.nanmean(np.array(LIST[ii+1])[maskwmo])
        std[ii, iw] = np.nanstd(np.array(LIST[ii+1])[maskwmo])

for ii in range(nL):
    pl.sca(axs[ii])
    pl.bar(LISTwmo, mean[ii,:], yerr=std[ii,:])
    pl.title('RMSD N3n layer ' + titleDICT[ii] + ' [mmol/m^3]')
    pl.grid()

pl.tight_layout()

pl.savefig(OUTDIR + 'meanRMSD_layer.png')


fig1b, axs = pl.subplots(nL, 1, sharex=True, sharey=True, figsize=[12, 8])

meanNorm = np.zeros((nL, len(LISTwmo)))
meanNorm[:,:] = np.nan
for iw, wmo in enumerate(LISTwmo):
    maskwmo = np.array(LIST[0]) == wmo
    for ii in range(nL):
        meanNorm[ii, iw] = np.nanmean(np.array(LIST[ii+1])[maskwmo]/np.array(LIST[4])[maskwmo])

for ii in range(nL):
    pl.sca(axs[ii])
    pl.bar(LISTwmo, meanNorm[ii,:])
    pl.title('RMSD norm 0-200 N3n layer ' + titleDICT[ii] + ' [mmol/m^3]')
    pl.grid()

pl.tight_layout()

pl.savefig(OUTDIR + 'meanRMSDnorm_layer.png')




thresholds = {
    0: [aa/10. for aa in range(20,4,-5)],
    1: [aa/10. for aa in range(30,14,-5)],
    2: [aa/10. for aa in range(35,19,-5)],
}

nexc = {}
Np = {}
for wmo in LISTwmo:
    maskwmo = np.array(LIST[0]) == wmo
    Np[wmo] = np.sum(maskwmo)
    nexc[wmo] = {}
    for ii in range(nL):
        nexc[wmo][ii] = np.zeros(len(thresholds[ii]))
        nexc[wmo][ii][:] = np.nan
        values = np.array(LIST[ii+1])[maskwmo]
        for iit,tt in enumerate(thresholds[ii][1:]):
            H_tt = thresholds[ii][iit]
            nexc[wmo][ii][iit+1] = np.nansum((values>tt) & (values<H_tt))
        nexc[wmo][ii][0] = np.nansum(values>thresholds[ii][0])



fig2, axs = pl.subplots(nL, 1, sharex=True, sharey=True, figsize=[12, 8])

for ii in range(nL):
    height = np.zeros((len(LISTwmo),len(thresholds[ii])))
    base = np.zeros((len(LISTwmo),len(thresholds[ii])))
    height[:,:] = np.nan
    base[:,:] = np.nan
    base[:,0] = 0.
    for iw,wmo in enumerate(LISTwmo):
        for iit in range(len(thresholds[ii])):
            height[iw,iit] = 100.*nexc[wmo][ii][iit]/Np[wmo]
        for itt in range(1,len(thresholds[ii])):
            base[iw,itt] = np.nansum(height[iw,:itt])

    category_colors = pl.get_cmap('RdYlGn')(np.linspace(0.15, 0.85, len(thresholds[ii])))
    pl.sca(axs[ii])
    for iit in range(len(thresholds[ii])):
        pl.bar(LISTwmo, height[:,iit], \
            bottom=base[:,iit], color=category_colors[iit], \
            label=np.str(thresholds[ii][iit]))
    # for iit in range(1,len(thresholds[ii])):
    #     pl.plot(LISTwmo, base[:,iit],'o')
    #     pl.plot(LISTwmo, base[:,iit],'-', color='grey')

    pl.legend(loc='best')
    pl.title('Perc excluded layer ' + titleDICT[ii])
    pl.grid()

pl.tight_layout()

pl.savefig(OUTDIR + 'percexclusion_layer.png')



thresholdsN = {
    0: [aa/10. for aa in range(40,9,-5)],
    1: [aa/10. for aa in range(50,19,-5)],
    2: [aa/10. for aa in range(50,19,-5)],
}

fixedThreshold = [3,4,5]

nexcN = {}
excprof = {}
for wmo in LISTwmo:
    maskwmo = np.array(LIST[0]) == wmo
    nexcN[wmo] = {}
    excprof[wmo] = np.zeros(Np[wmo])
    for ii in range(nL):
        nexcN[wmo][ii] = np.zeros(len(thresholdsN[ii]))
        nexcN[wmo][ii][:] = np.nan
        values = (np.array(LIST[ii+1])/np.array(LIST[4]))[maskwmo]
        for iit,tt in enumerate(thresholdsN[ii][1:]):
            H_tt = thresholdsN[ii][iit]
            nexcN[wmo][ii][iit+1] = np.nansum((values>tt) & (values<H_tt))
        nexcN[wmo][ii][0] = np.nansum(values>thresholdsN[ii][0])

        mask_fixt = values>fixedThreshold[ii]
        excprof[wmo][mask_fixt] += 2**ii 


fig2b, axs = pl.subplots(nL, 1, sharex=True, sharey=True, figsize=[12, 8])

for ii in range(nL):
    height = np.zeros((len(LISTwmo),len(thresholdsN[ii])))
    base = np.zeros((len(LISTwmo),len(thresholdsN[ii])))
    height[:,:] = np.nan
    base[:,:] = np.nan
    base[:,0] = 0.
    for iw,wmo in enumerate(LISTwmo):
        for iit in range(len(thresholdsN[ii])):
            height[iw,iit] = 100.*nexcN[wmo][ii][iit]/Np[wmo]
        for itt in range(1,len(thresholdsN[ii])):
            base[iw,itt] = np.nansum(height[iw,:itt])

    category_colors = pl.get_cmap('RdYlGn')(np.linspace(0.15, 0.85, len(thresholdsN[ii])))
    pl.sca(axs[ii])
    for iit in range(len(thresholdsN[ii])):
        pl.bar(LISTwmo, height[:,iit], \
            bottom=base[:,iit], color=category_colors[iit], \
            label=np.str(thresholdsN[ii][iit]))
    # for iit in range(1,len(thresholds[ii])):
    #     pl.plot(LISTwmo, base[:,iit],'o')
    #     pl.plot(LISTwmo, base[:,iit],'-', color='grey')

    pl.legend(loc='best')
    pl.title('Perc excluded layer ' + titleDICT[ii])
    pl.grid()

pl.tight_layout()

pl.savefig(OUTDIR + 'percexclusionNorm_layer.png')


percexc = np.zeros(len(LISTwmo))
percexc[:] = np.nan
for iw,wmo in enumerate(LISTwmo):
    percexc[iw] = 100.*np.nansum(excprof[wmo]>0)/Np[wmo]


wlayerDICT= {
    1: '200',
    2: '2-400',
    4: '4-600',
    3: '200 2-400',
    6: '2-400 4-600',
    5: '200 4-600',
    7: 'all',
}

layersH = np.zeros((len(LISTwmo),7))
layersH[:,:] = np.nan
layersB = np.zeros((len(LISTwmo),7))
layersB[:,:] = np.nan
layersB[:,0] = 0.
for iw,wmo in enumerate(LISTwmo):
    for ir in wlayerDICT.keys():
        if np.nansum(excprof[wmo]>0)==0: continue
        layersH[iw,ir-1] = 100.*np.nansum(excprof[wmo]==ir)/np.nansum(excprof[wmo]>0)
    for irr in range(1,len(wlayerDICT.keys())):
        layersB[iw,irr] = np.nansum(layersH[iw,:irr])


fig3, axs = pl.subplots(2, 1, sharex=True, sharey=True, figsize=[12, 8])


pl.sca(axs[0])
pl.bar(LISTwmo,percexc)
pl.title('Perc exclusion')
pl.grid()

pl.sca(axs[1])
for ir in range(len(wlayerDICT.keys())):
    pl.bar(LISTwmo,layersH[:,ir], \
        bottom=layersB[:,ir], \
        label=wlayerDICT[ir+1])

pl.legend(loc='best')
pl.title('Perc exclusion layer(s)')
pl.grid()


pl.tight_layout()

pl.savefig(OUTDIR + 'fixedt_perc.png')


pl.show(block=False)


