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

PATH=os.getcwd()   
SIM_NAME = PATH.strip('galileo/home/userexternal/eterzic0/BIOPTIMOD/.../bit.sea/')

maskfile    = '/galileo/home/userexternal/eterzic0/BIOPTIMOD/REA_16_T0/TEST01/wrkdir/bin/bit.sea.modified/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
jpk = nav_lev.shape[0]
mydepth = nav_lev[0:jpk-1]  
ncIN.close()


argslist = ['Ed380f', 'Ed412f', 'Ed490f']
kdlist   = ['Kd 380' , 'Kd 412' , 'Kd 490']

fig, ax = plt.subplots (nrows=3, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
fig.set_size_inches(12,12)

for i in range(len(argslist)):

    args=argslist[i]

    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
    TI = TimeInterval("20120101", "20171231","%Y%m%d")

    varname=var_conversions.FLOAT_OPT_BIOPTIMOD[args] 

    print('Calculating matchup for ', varname)

    MONTHS = np.arange(1, 12 + 1)
    months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

    FLOAT_W_mean = [] ; FLOAT_W_std = [] ; FLOAT_E_mean = [] ; FLOAT_E_std = []
    MODEL_W_mean = [] ; MODEL_W_std = [] ; MODEL_E_mean = [] ; MODEL_E_std = []

    for iMonth, month in enumerate(MONTHS):
        CLIM_MONTH_req       = timerequestors.Clim_month(month)

        print('MONTH:                            ', month)

        FLOAT_W, MODEL_W = calc_KD(varname, CLIM_MONTH_req, OGS.wes,  TI, M, mydepth, args)
        FLOAT_E, MODEL_E = calc_KD(varname, CLIM_MONTH_req, OGS.eas,  TI, M, mydepth, args)

        FLOAT_W_mean.append( np.mean(FLOAT_W)) ; FLOAT_E_mean.append( np.mean(FLOAT_E))
        FLOAT_W_std.append(  np.std( FLOAT_W)) ; FLOAT_E_std.append(  np.std( FLOAT_E))
        MODEL_W_mean.append( np.mean(MODEL_W)) ; MODEL_E_mean.append( np.mean(MODEL_E))
        MODEL_W_std.append(  np.std( MODEL_W)) ; MODEL_E_std.append(  np.std( MODEL_E))


    ax[i].scatter(MONTHS-0.25, FLOAT_W_mean, s=15, color='darkblue'    ,     label='float W')
    ax[i].scatter(MONTHS-0.15, MODEL_W_mean, s=15, color='dodgerblue'  ,     label='model W') 

    ax[i].scatter(MONTHS+0.15, FLOAT_E_mean, s=15, color='purple'     ,     label='float E')
    ax[i].scatter(MONTHS+0.25, MODEL_E_mean, s=15, color='palevioletred'  ,     label='model E')   

    ax[i].errorbar(MONTHS-0.25, FLOAT_W_mean, yerr=FLOAT_W_std, color='darkblue'   , fmt='o')
    ax[i].errorbar(MONTHS-0.15, MODEL_W_mean, yerr=MODEL_W_std, color='dodgerblue' , fmt='o')

    ax[i].errorbar(MONTHS+0.15, FLOAT_E_mean, yerr=FLOAT_E_std, color='purple'    , fmt='o')
    ax[i].errorbar(MONTHS+0.25, MODEL_E_mean, yerr=MODEL_E_std, color='palevioletred' , fmt='o')

    ax[i].set_xticks(MONTHS)
    ax[i].set_xticklabels(months_str)
    ax[i].set_ylabel('Kd [$m^{-1}$]')

    ax[i].tick_params(axis='both', which='major', labelsize=10)

    ax[0].legend(loc='upper center', ncol=2, fontsize=12)
    ax[i].set_title(kdlist[i], fontsize=16)

    countE, sigmaE, bias_valE, corr_coeffE, rE, pE, float_meanE, model_meanE = calc_statistics(np.array(FLOAT_E_mean), np.array(MODEL_E_mean))
    countW, sigmaW, bias_valW, corr_coeffW, rW, pW, float_meanW, model_meanW = calc_statistics(np.array(FLOAT_W_mean), np.array(MODEL_W_mean))
    
    file_dir = 'STATS/'
    
    file_out_E = file_dir +  args  + '_' + SIM_NAME.replace('/', '_') + 'E_monthly.stat'
    file_out_W = file_dir +  args  + '_' + SIM_NAME.replace('/', '_') + 'W_monthly.stat'
    
    f_out_E   = writefile(file_out_E, args, countE, sigmaE, corr_coeffE, bias_valE, rE, pE, float_meanE, model_meanE)
    f_out_W   = writefile(file_out_W, args, countW, sigmaW, corr_coeffW, bias_valW, rW, pW, float_meanW, model_meanW)

fig.suptitle(SIM_NAME, fontsize=20)
plot_out = 'PLOTS/' + SIM_NAME.replace('/', '_') + '_scatter_MONTHLY_WE.png'
fig.savefig(plot_out, format='png',dpi=150)

print('Calculation successfully computed.')
