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
SIM_NAME = os.path.join(*PATH.split('/')[-3:-1])#PATH.strip('galileo/home/userexternal/eterzic0/BIOPTIMOD/DCM_VAL/INPUT/.../bit.sea/')

maskfile    = '/galileo/home/userexternal/eterzic0/BIOPTIMOD/KD_VAL/INPUT/REA_16_T0/TEST01/wrkdir/bin/bit.sea.modified/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
jpk = nav_lev.shape[0]
mydepth = nav_lev[0:jpk-1]  
ncIN.close()

varlist   = ['DCM depth' , 'Chl at DCM' ]
arglist   = ['DCM', 'CHL']
labellist = ['DCM depth [m]', 'Chl at DCM [$mg \, m^{-3} $]']

for i in range(len(varlist)):

    fig, ax = plt.subplots (nrows=1, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
    fig.set_size_inches(12,4)

    args='P_l'

    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
    TI = TimeInterval("20120101", "20171231","%Y%m%d")

    varname=var_conversions.FLOAT_OPT_BIOPTIMOD[args] 

    print('Calculating matchup for ', varname)

    MONTHS = np.arange(1, 12 + 1)
    months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

    FLOAT_W_median = [] ; FLOAT_W_IQR = [] ; FLOAT_E_median = [] ; FLOAT_E_IQR = []
    MODEL_W_median = [] ; MODEL_W_IQR = [] ; MODEL_E_median = [] ; MODEL_E_IQR = []

    FLOAT_W_25 = [] ; FLOAT_E_25 = [] ; FLOAT_W_75 = [] ; FLOAT_E_75 = [] ; 
    MODEL_W_25 = [] ; MODEL_E_25 = [] ; MODEL_W_75 = [] ; MODEL_E_75 = []

    FLOAT_W_5 = [] ; FLOAT_E_5 = [] ; FLOAT_W_95 = [] ; FLOAT_E_95 = [] ; 
    MODEL_W_5 = [] ; MODEL_E_5 = [] ; MODEL_W_95 = [] ; MODEL_E_95 = []

    for iMonth, month in enumerate(MONTHS):
        CLIM_MONTH_req       = timerequestors.Clim_month(month)

        print('MONTH:                            ', month)

        DCM_FLOAT_W, DCM_MODEL_W, CHL_FLOAT_W, CHL_MODEL_W = calc_DCM(varname, CLIM_MONTH_req, OGS.wes,  TI, M, mydepth, args)
        DCM_FLOAT_E, DCM_MODEL_E, CHL_FLOAT_E, CHL_MODEL_E = calc_DCM(varname, CLIM_MONTH_req, OGS.eas,  TI, M, mydepth, args)

        if i == 0: 
            FLOAT_W = DCM_FLOAT_W ; MODEL_W = DCM_MODEL_W ; FLOAT_E = DCM_FLOAT_E ; MODEL_E = DCM_MODEL_E
        else:
            FLOAT_W = CHL_FLOAT_W ; MODEL_W = CHL_MODEL_W ; FLOAT_E = CHL_FLOAT_E ; MODEL_E = CHL_MODEL_E

        FLOAT_W_median.append( np.nanmedian(FLOAT_W)) ; FLOAT_E_median.append( np.nanmedian(FLOAT_E))

        # All the percentiles
        FLOAT_W_5.append(np.nanpercentile( FLOAT_W,  5)) ; FLOAT_E_5.append(np.nanpercentile( FLOAT_E,  5))
        FLOAT_W_25.append(np.nanpercentile(FLOAT_W, 25)) ; FLOAT_E_25.append(np.nanpercentile(FLOAT_E, 25))
        FLOAT_W_75.append(np.nanpercentile(FLOAT_W, 75)) ; FLOAT_E_75.append(np.nanpercentile(FLOAT_E, 75))
        FLOAT_W_95.append(np.nanpercentile(FLOAT_W, 95)) ; FLOAT_E_95.append(np.nanpercentile(FLOAT_E, 95))

        MODEL_W_median.append( np.nanmedian(MODEL_W)) ; MODEL_E_median.append( np.nanmedian(MODEL_E))

        # All the percentiles
        MODEL_W_5.append(np.nanpercentile( MODEL_W,  5)) ; MODEL_E_5.append(np.nanpercentile( MODEL_E,  5))
        MODEL_W_25.append(np.nanpercentile(MODEL_W, 25)) ; MODEL_E_25.append(np.nanpercentile(MODEL_E, 25))
        MODEL_W_75.append(np.nanpercentile(MODEL_W, 75)) ; MODEL_E_75.append(np.nanpercentile(MODEL_E, 75))
        MODEL_W_95.append(np.nanpercentile(MODEL_W, 95)) ; MODEL_E_95.append(np.nanpercentile(MODEL_E, 95))

    FLOAT_W_IQR = 0.5*(np.array(FLOAT_W_75) - np.array(FLOAT_W_25))    
    FLOAT_E_IQR = 0.5*(np.array(FLOAT_E_75) - np.array(FLOAT_E_25))
    MODEL_W_IQR = 0.5*(np.array(MODEL_W_75) - np.array(MODEL_W_25))    
    MODEL_E_IQR = 0.5*(np.array(MODEL_E_75) - np.array(MODEL_E_25))

    ax.scatter(MONTHS-0.25, FLOAT_W_median, s=15, color='darkblue'    ,     label='float W')
    ax.scatter(MONTHS-0.15, MODEL_W_median, s=15, color='dodgerblue'  ,     label='model W') 

    ax.scatter(MONTHS+0.15, FLOAT_E_median, s=15, color='purple'         ,     label='float E')
    ax.scatter(MONTHS+0.25, MODEL_E_median, s=15, color='palevioletred'  ,     label='model E')   

    ax.errorbar(MONTHS-0.25, FLOAT_W_median, yerr=FLOAT_W_IQR, color='darkblue'   , fmt='o')
    ax.errorbar(MONTHS-0.15, MODEL_W_median, yerr=MODEL_W_IQR, color='dodgerblue' , fmt='o')

    ax.errorbar(MONTHS+0.15, FLOAT_E_median, yerr=FLOAT_E_IQR, color='purple'        , fmt='o')
    ax.errorbar(MONTHS+0.25, MODEL_E_median, yerr=MODEL_E_IQR, color='palevioletred' , fmt='o')

    ax.scatter(MONTHS-0.25, FLOAT_W_5, s=10, color='darkblue'  )
    ax.scatter(MONTHS-0.25, FLOAT_W_95, s=10, color='darkblue'  )

    ax.scatter(MONTHS-0.15, MODEL_W_5, s=10, color='dodgerblue'  )
    ax.scatter(MONTHS-0.15, MODEL_W_95, s=10, color='dodgerblue'  )

    ax.scatter(MONTHS+0.15, FLOAT_E_5, s=10, color='purple'  )
    ax.scatter(MONTHS+0.15, FLOAT_E_95, s=10, color='purple'  )

    ax.scatter(MONTHS+0.25, MODEL_E_5, s=10, color='palevioletred'  )  
    ax.scatter(MONTHS+0.25, MODEL_E_95, s=10, color='palevioletred'  )

    ax.set_xticks(MONTHS)
    ax.set_xticklabels(months_str)
    ax.set_ylabel(labellist[i])

    ax.tick_params(axis='both', which='major', labelsize=10)

    ax.legend(loc='upper right', ncol=2, fontsize=12)
    ax.set_title(varlist[i], fontsize=16)

    countE, sigmaE, bias_valE, corr_coeffE, rE, pE, float_meanE, model_meanE = calc_statistics(np.array(MODEL_E_median), np.array(FLOAT_E_median))
    countW, sigmaW, bias_valW, corr_coeffW, rW, pW, float_meanW, model_meanW = calc_statistics(np.array(MODEL_W_median), np.array(FLOAT_W_median))
    
    file_dir = 'STATS_MEDIAN/'
    
    file_out_E = file_dir +  arglist[i]  + '_' + SIM_NAME.replace('/', '_') + 'E_monthly.stat'
    file_out_W = file_dir +  arglist[i]  + '_' + SIM_NAME.replace('/', '_') + 'W_monthly.stat'
    
    f_out_E   = writefile(file_out_E, arglist[i], countE, sigmaE, corr_coeffE, bias_valE, rE, pE, float_meanE, model_meanE)
    f_out_W   = writefile(file_out_W, arglist[i], countW, sigmaW, corr_coeffW, bias_valW, rW, pW, float_meanW, model_meanW)

    #fig.suptitle(SIM_NAME, fontsize=20)
    plot_out = 'PLOTS_MEDIAN/' + SIM_NAME.replace('/', '_') + '_' + arglist[i] + '_scatter_MONTHLY_WE.png'
    fig.savefig(plot_out, format='png',dpi=150)

print('Calculation successfully computed.')
