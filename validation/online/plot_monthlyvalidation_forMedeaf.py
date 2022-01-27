# Plot results of Monthly Validation Report 
# taken from MVR file sent to CMEMS.
# For sat and float

import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    plot something
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output dir'''
                            )

    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required = True,
                            default = '',
                            help = 'Input dir')


    return parser.parse_args()
args = argument()


import numpy as np
import netCDF4 as NC
import datetime
import matplotlib.pyplot as plt

from commons.utils import addsep
from commons.Timelist import TimeList
from basins import V2 as OGS



INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)


TLmvr = TimeList.fromfilenames(None,INDIR,'product_quality*nc',
                        prefix='product_quality_stats_MEDSEA_ANALYSIS_FORECAST_BIO_006_014_',
                        dateformat='%Y%m%d')

dates = []
satstats = []
floatstats = []

DICTvardim = {
    'areas': 'area_names',
    'metrics': 'metric_names',
    'forecasts': 'forecasts',
    'time': 'time',
    'depths': 'depths',
}

for ii,filein in enumerate(TLmvr.filelist):
    print filein
    MVR = NC.Dataset(filein,'r')
    datesmonth = MVR.variables['time'][:].data.copy()
    dates.extend(list(datesmonth))
    satstats_month = MVR.variables['stats_surface_chlorophyll'][:].data.copy()
    satstats.extend(list(satstats_month))
    floatstats_month = MVR.variables['stats_profile_chlorophyll'][:].data.copy()
    floatstats.extend(list(floatstats_month))
    if ii==0:
        DICTdim_sat = {}
        dimlist = [dd.encode() for dd in MVR.variables['stats_surface_chlorophyll'].dimensions]
        for iid,dd in enumerate(dimlist):
            if 'surf' in dd:
                DICTdim_sat[dd] = ['surface',iid]
                continue
            varname = DICTvardim[dd]
            vv = MVR.variables[varname][:].data.copy()
            if 'float' in vv.dtype.name: 
                DICTdim_sat[dd] = [vv,iid]
            if 'string' in vv.dtype.name: 
                vLIST = []
                for iiv in range(vv.shape[0]):
                    vLIST.append(''.join(vv[iiv,:]))
                DICTdim_sat[dd] = [vLIST,iid] 
        
        DICTdim_float = {}
        dimlist = [dd.encode() for dd in MVR.variables['stats_profile_chlorophyll'].dimensions]
        for iid,dd in enumerate(dimlist):
            varname = DICTvardim[dd]
            vv = MVR.variables[varname][:].data.copy()
            if 'float' in vv.dtype.name: 
                DICTdim_float[dd] = [vv,iid]
            if 'string' in vv.dtype.name: 
                vLIST = []
                for iiv in range(vv.shape[0]):
                    vLIST.append(''.join(vv[iiv,:]))
                DICTdim_float[dd] = [vLIST,iid] 
    # break


array_floatstats = np.array(floatstats)
array_floatstats[array_floatstats>1.e+19] = np.nan
array_satstats = np.array(satstats)
array_satstats[array_satstats>1.e+19] = np.nan

dates_datetime = []
for dd in dates:
    ddordinal = np.int(dd) + datetime.datetime(1970,1,1).toordinal()
    dd_datetime = datetime.datetime.fromordinal(ddordinal)
    dates_datetime.append(dd_datetime)


DICTsubgroup_index = {}
for subaggregate in OGS.MVR.basin_list:
    for isub,sub in enumerate(OGS.P.basin_list):
        if subaggregate.name in sub.name:
            DICTsubgroup_index[subaggregate.name] = isub
            break

DICTsub_shortname = {}
for subname in DICTdim_sat['areas'][0]:
    for sub in OGS.P.basin_list:
        if subname in sub.extended_name:
            DICTsub_shortname[subname] = sub.name
            break

DICTalpha = {
    12: 1,
    36: 0.6,
    60: 0.3,
    -12: 1,
}

DICTvargroup = {
    'number of data values': 0,
    'mean of product':       1,
    'mean of reference':     1,
    'mean squared error':    2,
    'variance of product':   3,
    'variance of reference': 3,
    'correlation':           4,
    'anomaly correlation':   5,
}


cmap = plt.get_cmap("Dark2")

plt.close('all')
indmetrics = DICTdim_sat['metrics'][1]
indsub = DICTdim_sat['areas'][1]

noforecasts = ['number of data values','mean of reference','variance of reference']
for isub,subname in enumerate(DICTdim_sat['areas'][0]):
    print (subname)
    print ('...........SAT.....')
    fig,axs = plt.subplots(3,2,sharex=True,figsize=[14,8])#,sharey=True)
    for iim,mm in enumerate(DICTdim_sat['metrics'][0]):
        print (mm)
        nax = DICTvargroup[mm]
        ix_ax = nax/2
        iy_ax = nax-2*ix_ax
        plt.sca(axs[ix_ax,iy_ax])
        for iif,ff in enumerate(DICTdim_sat['forecasts'][0]):
            if 'reference' in mm:
                label = mm
            else:
                label = mm + ' ' + np.str(ff)
            if ff<0:
                linestyle=':'
            else:
                linestyle='-'
            
            DICTind = {
                indmetrics: iim,
                indsub: isub,
                DICTdim_sat['forecasts'][1]: iif,
            }
            selection = [DICTind.get(dd,slice(None)) for dd in range(len(DICTdim_sat.keys()))]
            plt.plot(dates_datetime,array_satstats[selection],
                            color=cmap(iim),
                            linestyle=linestyle,
                            alpha=DICTalpha[ff],
                            label=label)
            if mm in noforecasts: break


    for axx in axs.flatten():    
        plt.sca(axx)
        plt.legend()
        plt.grid()
    plt.suptitle(subname)
    fig.autofmt_xdate()

    # plt.show(block=False)
    plt.savefig(OUTDIR + '/satmetric' + DICTsub_shortname[subname].upper() + '.png')




# float
plt.close('all')
indmetrics = DICTdim_sat['metrics'][1]
indsub = DICTdim_sat['areas'][1]

LISTalpha_depth = [1-x*2/10. for x in range(len(DICTdim_float['depths'][0]))]

for sub in OGS.MVR.basin_list:
    print (sub.name)
    print ('..........FLOAT......')
    fig,axs = plt.subplots(3,2,sharex=True,figsize=[14,8])#,sharey=True)
    for iim,mm in enumerate(DICTdim_sat['metrics'][0]):
        print (mm)
        nax = DICTvargroup[mm]
        ix_ax = nax/2
        iy_ax = nax-2*ix_ax
        plt.sca(axs[ix_ax,iy_ax])
        label = mm
        for iid,depth in enumerate(DICTdim_float['depths'][0]):
            if 'reference' in mm:
                label = mm
            if iid>0:
                label = None
            DICTind = {
                indmetrics: 0,
                indsub: DICTsubgroup_index[sub.name],
                DICTdim_sat['forecasts'][1]: 0,
                DICTdim_float['depths'][1]: iid,
            }
            selection = [DICTind.get(dd,slice(None)) for dd in range(len(DICTdim_float.keys()))]
            # if np.nansum(array_floatstats[selection])==0: continue


            DICTind = {
                indmetrics: iim,
                indsub: DICTsubgroup_index[sub.name],
                DICTdim_sat['forecasts'][1]: 0,
                DICTdim_float['depths'][1]: iid,
            }
            selection = [DICTind.get(dd,slice(None)) for dd in range(len(DICTdim_float.keys()))]
            plt.plot(dates_datetime,array_floatstats[selection],
                            color=cmap(iim),
                            alpha=LISTalpha_depth[iid],
                            label=label)


    for axx in axs.flatten():
        plt.sca(axx)
        plt.legend()
        plt.grid()
    plt.suptitle(sub.name)
    fig.autofmt_xdate()

    # plt.show(block=False)
    plt.savefig(OUTDIR + '/floatmetric' + sub.name.upper() + '.png')
