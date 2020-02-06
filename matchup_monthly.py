from ancillary import *
from basins import V2 as OGS
from commons.layer import Layer
from commons import timerequestors
from instruments import optbio_float_2019
from instruments import var_conversions
import matplotlib.pyplot as plt
import numpy as np
import os,sys
from profiler_ET import *
import scipy.io.netcdf as NC



maskfile    = '/galileo/home/userexternal/eterzic0/BIOPTIMOD/REA_16_T0/TEST01/wrkdir/bin/bit.sea.modified/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
jpk = nav_lev.shape[0]
mydepth = nav_lev[0:jpk-1]  
ncIN.close()


argslist = ['Ed380f', 'Ed412f', 'Ed490f']
kdlist   = ['Kd 380' , 'Kd 412' , 'Kd 490']

fig, ax = plt.subplots (nrows=3, ncols=1)
fig.set_size_inches(8,12)

for i in range(len(argslist)):

    args=argslist[i]

    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
    TI = TimeInterval("20120101", "20171231","%Y%m%d")

    varname=var_conversions.FLOAT_OPT_BIOPTIMOD[args] 

    print('Calculating matchup for ', varname)

    MONTHS = np.arange(1, 12 + 1)
    months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

    FLOAT_mean = [] ; FLOAT_std = []
    MODEL_mean = [] ; MODEL_std = []

    for iMonth, month in enumerate(MONTHS):
        CLIM_MONTH_req       = timerequestors.Clim_month(month)

        print('MONTH:                            ', month)

        FLOAT, MODEL = calc_KD(varname, CLIM_MONTH_req, OGS.Pred,  TI, M, mydepth, args)

        FLOAT_mean.append( np.mean(FLOAT)) 
        FLOAT_std.append(  np.std( FLOAT))

        MODEL_mean.append( np.mean(MODEL)) 
        MODEL_std.append(  np.std( MODEL))


    ax[i].scatter(MONTHS-0.1, FLOAT_mean, s=100, color='darkblue', label='float')
    ax[i].scatter(MONTHS+0.1, MODEL_mean, s=100, color='dodgerblue', label='model')   

    ax[i].errorbar(MONTHS-0.1, FLOAT_mean, yerr=FLOAT_std, color='darkblue', fmt='o')
    ax[i].errorbar(MONTHS+0.1, MODEL_mean, yerr=MODEL_std, color='dodgerblue', fmt='o')

    ax[i].set_xticks(MONTHS)
    ax[i].set_xticklabels(months_str)

    ax[i].tick_params(axis='both', which='major', labelsize=10)

    ax[0].legend(loc='upper center', fontsize=16)
    ax[i].set_title(kdlist[i], fontsize=18)

    count, sigma, bias_val, corr_coeff, r, p, float_mean, model_mean = calc_statistics(np.array(FLOAT_mean), np.array(MODEL_mean)

    file_dir = PATH + 'STATS/'
    file_out =  file_dir +  args  + '_monthly.stat'
    
    f_out   = writefile(file_out, args, count, sigma, corr_coeff, bias_val, r, p, float_mean, model_mean)


plot_out = 'PLOTS/scatter_MONTHLY.png'
fig.savefig(plot_out, format='png',dpi=150)

print('Calculation successfully computed.')
