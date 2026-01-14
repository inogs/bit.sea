import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot results of Monthly Validation Report
    taken from MVR file sent to CMEMS.
    For sat and float
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

    parser.add_argument(   '--varname', '-vv',
                            type = str,
                            required = True,
                            default = 'chlorophyll',
                            help = 'varname')


    return parser.parse_args()
args = argument()


import numpy as np
import netCDF4 as NC
import datetime
import matplotlib.pyplot as plt

from bitsea.commons.utils import addsep
from bitsea.commons.Timelist import TimeList
from bitsea.basins import V2 as OGS
import sys
import pandas as pd
import os

INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)


TLmvr = TimeList.fromfilenames(None,INDIR,'product_quality*nc',
                        prefix='product_quality_stats_MEDSEA_ANALYSISFORECAST_BGC_006_014_',
                        dateformat='%Y%m%d')

dates = []
satstats = []
floatstats = []

VAR = args.varname
LIST_VAR = ['chlorophyll', 'nitrate', 'oxygen']

if VAR not in LIST_VAR:
    sys.exit(f"Selected varname '{VAR}' is not valid. Choose from: {', '.join(LIST_VAR)}")


DICTvardim = {
    'time': 'time',    
    'area': 'area',
    'metric': 'metric',
    'depth': 'depth',
    'layer': 'layer',
    'forecast': 'forecast',
}


def reshape_label(handles_labels):
   """ Dividere la lista in due sotto-liste: una per "Mod" e una per "Ref
   """
   mod_list = [item for item in handles_labels if 'Mod' in item[0]]
   ref_list = [item for item in handles_labels if 'Ref' in item[0]]
   reshaped_list = []
   for mod, ref in zip(mod_list, ref_list):
       reshaped_list.append(mod)
       reshaped_list.append(ref)
   return(reshaped_list)


def plot_profiles_rmse(df, outdir, VAR, sub):
    """Crea e salva un line plot con doppio asse Y:
    """
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax2 = ax1.twinx()

    # Linee
    ax1.plot(df.index, df['#profiles'],
             linewidth=3, marker='o',
             color='tab:blue', label='#profiles')

    ax2.plot(df.index, df['MSE'],
             linewidth=3, marker='s',
             color='tab:orange', label='MSE')

    # Assi
    ax1.set_xlabel('Depth [m]', fontsize=12)
    ax1.set_ylabel('#profiles', fontsize=12, color='tab:blue')
    ax2.set_ylabel('MSE', fontsize=12, color='tab:orange')

    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax2.tick_params(axis='y', labelcolor='tab:orange')

    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'{VAR} – {sub.name.upper()}', fontsize=14)

    # Legenda unica
    l1, lab1 = ax1.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lab1 + lab2, loc='best')

    plt.tight_layout()

    outfile = os.path.join(
        outdir,f"{VAR}_floatmetric_{sub.name.upper()}.png")

    plt.savefig(outfile)
    plt.close(fig)


for ii,filein in enumerate(TLmvr.filelist):
    print(filein)
    MVR = NC.Dataset(filein,'r')
    datesmonth = MVR.variables['time'][:].data.copy()
    dates.extend(list(datesmonth))
    if VAR == "chlorophyll":
       satstats_month = MVR.variables['stats_chlorophyll-a_sat-l3'][:].data.copy()
       satstats.extend(list(satstats_month))
       print('stats_'+VAR+'-a-ins-pf')
       floatstats_month = MVR.variables['stats_'+VAR+'-a_ins-pf'][:].data.copy()
    #o2;no3# 
    else:
       NAMEVAR='stats_'+VAR+'-ins-pf'
       print(NAMEVAR)
       floatstats_month = MVR.variables['stats_'+VAR+'_ins-pf'][:].data.copy()
    floatstats.extend(list(floatstats_month))
    if ii==0:
        DICTdim_sat = {}
        if VAR == "chlorophyll":
           dimtuple = MVR.variables['stats_chlorophyll-a_sat-l3'].dimensions
           UNITS_for_Title=("(mg Chl*m$^{-3}$)")
        else: 
           dimtuple = MVR.variables['stats_'+VAR+'_ins-pf' ].dimensions 
           if VAR=="nitrate": UNITS_for_Title = '(mmol NO3*m$^{-3}$)'
           else: UNITS_for_Title = '(mmol O2*m$^{-3})$'

        for iid,dd in enumerate(dimtuple):
           if VAR == "chlorophyll":
              if 'depth' in dd:
                 DICTdim_sat[dd] = ['depth',iid]
                 continue
           else:
              if 'lay' in dd: 
                 DICTdim_sat[dd] = ['layer',iid]
                 continue 
           varname = DICTvardim[dd]
           vv = MVR.variables[varname][:]
           if vv.dtype.kind=='f': 
               DICTdim_sat[dd] = [vv,iid]
           if vv.dtype.kind=='S' or vv.dtype.kind=='O': 
               vLIST = []
               #for iiv in range(vv.shape[0]):
                   #vLIST.append(''.join([vv[iiv,kk].decode("utf-8") for kk in range(vv.shape[1])]))
               vLIST=list(vv)
               DICTdim_sat[dd] = [vLIST,iid] 
       
        DICTdim_float = {}
        if VAR == "chlorophyll":
            dimtuple = MVR.variables['stats_'+VAR+'-a_ins-pf' ].dimensions
        else:
            dimtuple = MVR.variables['stats_'+VAR+'_ins-pf' ].dimensions
        for iid,dd in enumerate(dimtuple):
            varname = DICTvardim[dd]
            vv = MVR.variables[varname][:]
            if vv.dtype.kind=='f': 
                DICTdim_float[dd] = [vv,iid]
            if vv.dtype.kind=='S' or vv.dtype.kind=='O': 
                vLIST = list(vv)
                #for iiv in range(vv.shape[0]):
                #    vLIST.append(''.join([vv[iiv,kk].decode("utf-8") for kk in range(vv.shape[1])]))
                DICTdim_float[dd] = [vLIST,iid] 


array_floatstats = np.array(floatstats)
array_floatstats[array_floatstats>1.e+19] = np.nan
array_satstats = np.array(satstats)
array_satstats[array_satstats>1.e+19] = np.nan

dates_datetime = []
for dd in dates:
    ddordinal = int(dd) + datetime.datetime(1950,1,1).toordinal()
    dd_datetime = datetime.datetime.fromordinal(ddordinal)
    dates_datetime.append(dd_datetime)


DICTsubgroup_index = {}
for subaggregate in OGS.MVR.basin_list:
    for isub,sub in enumerate(OGS.P.basin_list):
        if subaggregate.name in sub.name:
            DICTsubgroup_index[subaggregate.name] = isub
            break

DICTsub_shortname = {}
for subname in DICTdim_sat['area'][0]:
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
    'anomaly correlation':   4,
    'correlation':           5,
}


cmap = plt.get_cmap("Dark2")
cmap_edge = plt.get_cmap("gray")

plt.close('all')

if VAR == 'chlorophyll':

    indmetrics = DICTdim_sat['metric'][1]
    indsub = DICTdim_sat['area'][1]

    noforecasts = ['number of data values','mean of reference','variance of reference']
    for isub,subname in enumerate(DICTdim_sat['area'][0]):
       fig,axs = plt.subplots(3,2,sharex=True,figsize=[14,8])#,sharey=True)
       for iim,mm in enumerate(DICTdim_sat['metric'][0]):
           nax = DICTvargroup[mm]
           ix_ax = int(np.floor(nax/2))
           iy_ax = nax-2*ix_ax
           plt.sca(axs[ix_ax,iy_ax])
           for iif,ff in enumerate(DICTdim_sat['forecast'][0]):
               if 'reference' in mm:
                   label = mm
               else:
                   label = mm + ' ' + str(ff)
               if ff<0:
                   linestyle=':'
               else:
                   linestyle='-'
            
               DICTind = {
                   indmetrics: iim,
                   indsub: isub,
                   DICTdim_sat['forecast'][1]: iif,
               }
               selection_sat = [DICTind.get(dd,slice(None)) for dd in range(len(DICTdim_sat.keys()))]
               plt.plot(dates_datetime,array_satstats[tuple(selection_sat)],
                               color=cmap(iim),
                               linestyle=linestyle,
                               alpha=DICTalpha[ff],
                               label=label)
               if mm in noforecasts: break

       for axx in axs.flatten():    
           plt.sca(axx)
           plt.legend()
           plt.grid()
       plt.suptitle(subname + ' ' + UNITS_for_Title )
       fig.autofmt_xdate()
       BASIN=subname.replace(' ','_')
       plt.savefig(OUTDIR + '/satmetric_' +  BASIN+ '.png')

else:
    pass


# float
plt.close('all')
indmetrics = DICTdim_sat['metric'][1]
indsub     = DICTdim_sat['area'][1]
indlayer   = DICTdim_float['layer'][1]

LISTalpha_depth = np.linspace(1,0.3, len(DICTdim_float['layer'][0]))

for sub in OGS.MVR.basin_list:
    print ('Processing ' +  sub.name)
    df=pd.DataFrame(index=np.arange(0,7), columns=['#profiles','MSE'])
    for iid,depth in enumerate(DICTdim_float['layer'][0]):
        handles_labels = []
        skip_figure = False
        fig,axs = plt.subplots(3,2,sharex=True,figsize=[14,12])#,sharey=True)
        for iim,mm in enumerate(DICTdim_sat['metric'][0]): # loop su tutte le metriche 
            nax = DICTvargroup[mm]
            ix_ax = int(np.floor(nax/2))
            iy_ax = nax-2*ix_ax
            ax = axs[ix_ax, iy_ax]
            group_label = "Mod" if iim == 0 else "Ref"
            label = f"{depth}m {group_label}"
            lab   = str(int(depth))
            CMAP=cmap(0)
            title= mm.capitalize()
                
            if str(nax) in ['1','3']: # see DICTvargroup
                title= mm.split()[0].capitalize() + ' Mod vs Ref'
                if 'ref' in mm:
                    CMAP=cmap(1)

            #ax.set_title(title)

            DICTind = {
                indmetrics: iim,
                indsub: DICTsubgroup_index[sub.name],
                DICTdim_sat['forecast'][1]: 0,
                DICTdim_float['layer'][1]: iid,
            }
            selection = [DICTind.get(dd,slice(None)) for dd in range(len(DICTdim_float.keys()))]
            
            arr=array_floatstats[tuple(selection)]
            if mm.startswith('number'):
                if np.all(arr == 0) or np.all(np.isnan(arr)):
                   print('No data found in --> ' + sub.name + ' in depth -->' + str(int(depth)  ))
                   skip_figure = True
                   break 

                else:
                   arr[arr == 0] = np.nan
            
            MEAN_VAL= np.nanmean(arr)
            if mm.startswith('number') or mm == ('mean squared error'):
                ax.set_title(title + ' ( mean value: ' +str(MEAN_VAL)  +') ' )
            else:
                ax.set_title(title)

            if mm.startswith('number'):
               df.at[iid, '#profiles'] = np.around(MEAN_VAL, 2)
            elif mm == 'mean squared error':
               df.at[iid, 'MSE'] = MEAN_VAL

            # plot line and scatter dots
            ax.plot(dates_datetime,arr,color='k',linewidth=1,zorder=1) #alpha=LISTalpha_depth[iid])
            line = ax.scatter(dates_datetime,arr,facecolor=CMAP,    
                   edgecolors='k',label=label,s=70, marker='o',zorder=2)
            
            y_values = array_floatstats[tuple(selection)]
            x_values = dates_datetime

            if iim == 0 or iim==2:
               handles_labels.append((label, line))
        
        if skip_figure:
            plt.close(fig)
            continue

        for axx in axs.flatten():
            axx.grid(True, linestyle='--', linewidth=0.5)
        
        plt.suptitle(f"{sub.name.capitalize()} {VAR} {UNITS_for_Title} layer:  {int(depth)}m" )

        fig.autofmt_xdate()

        handles_labels=reshape_label(handles_labels)
        handles_labels_mod = []
        handles_labels_ref = []


        labels_seen = set()
        handles, labels = [], []
        for label, handle in handles_labels:
            if label not in labels_seen:
                labels_seen.add(label)
                labels.append(label)
                handles.append(handle)

        fig.legend(handles, labels,
                   loc='lower center',
                   ncol=2,
                   fontsize=18,
                   bbox_to_anchor=(0.5, -0.))  # Spazio sotto

        plt.subplots_adjust(top=0.93, left=0.04, bottom=0.12, right=0.97)
        plt.savefig(OUTDIR + '/'+VAR +'_floatmetric_' + sub.name.upper() +'_'+ str(int(depth)) + 'm.png')
        plt.close(fig)
    
    df.index= list(DICTdim_float['layer'][0])
    outdir = os.path.join(OUTDIR, 'CSV')
    os.makedirs(outdir, exist_ok=True)
    plot_profiles_rmse(df, outdir, VAR, sub)
    df.to_csv(outdir +'/'+ VAR +'_floatmetric_' + sub.name.upper() + '.csv' )
