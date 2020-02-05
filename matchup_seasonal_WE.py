from ancillary import *
import argparse
from basins import V2 as OGS
from commons.layer import Layer
from commons import timerequestors
from commons import season as season
from instruments import optbio_float_2019
from instruments import var_conversions
import matplotlib.pyplot as plt
import numpy as np
import os,sys
from profiler_ET import *
import scipy.io.netcdf as NC

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates scatter plots of surface irradiance data 
    from BGC-Argo floats and the atmospheric OASIM output
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--var', '-v',
                                type = str,
                                required = True,
                                choices = ['Ed380f','Ed412f','Ed490f'],
                                help = ''' Choose the variable to plot.'''
                                )

    return parser.parse_args()


args = argument()
#args = 'Ed490f'  
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

maskfile    = '/galileo/home/userexternal/eterzic0/BIOPTIMOD/REA_16_T0/TEST01/wrkdir/bin/bit.sea.modified/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
jpk = nav_lev.shape[0]
mydepth = nav_lev[0:jpk-1]  
ncIN.close()

TI = TimeInterval("20120101", "20171231","%Y%m%d")

#varname=var_conversions.FLOAT_OPT_BIOPTIMOD['Ed490f'] 
varname=var_conversions.FLOAT_OPT_BIOPTIMOD[args.var]

print('Calculating matchup for ', varname)

seasonobj = season.season()

SEAS_STR = ['WINTER', 'SPRING', 'SUMMER', 'AUTUMN']

fig,axs = plt.subplots(2,2, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
fig.set_size_inches(12,9)
PATH='/galileo/home/userexternal/eterzic0/BIOPTIMOD/REA_16_INIT/TEST03/bit.sea/'
SIM_NAME = PATH.strip('galileo/home/userexternal/eterzic0/BIOPTIMOD/.../bit.sea/')
for iseas in range(len(SEAS_STR)):  
    print('SEASON: ', SEAS_STR[iseas])
    seasonobj = season.season()
    CLIM_SEAS_req = timerequestors.Clim_season(iseas, seasonobj)

    FLOAT_W, MODEL_W = calc_KD(varname, CLIM_SEAS_req, OGS.wes,  TI, M, mydepth, args.var)
    FLOAT_E, MODEL_E = calc_KD(varname, CLIM_SEAS_req, OGS.eas,  TI, M, mydepth, args.var)

    count1, sigma1, corr_coeff1, bias_val1, b1, a1, float_mean1, model_mean1 = calc_statistics(FLOAT_W, MODEL_W) 
    count2, sigma2, corr_coeff2, bias_val2, b2, a2, float_mean2, model_mean2 = calc_statistics(FLOAT_E, MODEL_E)

    ax1 = plot_basin(iseas, axs, FLOAT_W, MODEL_W, 'darkblue',   OGS.wes, SEAS_STR, 0.05, 0.95,sigma1, bias_val1, corr_coeff1, b1, a1, count1)
    ax2 = plot_basin(iseas, axs, FLOAT_E, MODEL_E, 'dodgerblue', OGS.eas, SEAS_STR, 0.75, 0.35,sigma2, bias_val2, corr_coeff2, b2, a2, count2)
    
    plot_out = PATH + 'PLOTS/scatter_' + args.var + '_SEAS_WE.png'
    fig.suptitle(SIM_NAME + ' ' + varname)
    fig.savefig(plot_out, format='png',dpi=150)

    file_dir = PATH + 'STATS/'
    file_out1 =  file_dir +  args.var + '_' + SEAS_STR[iseas] + '_WEST.stat'
    file_out2 =  file_dir +  args.var + '_' + SEAS_STR[iseas] + '_EAST.stat'
    f_out1    = writefile(file_out1, args.var, count1, sigma1, corr_coeff1, bias_val1, b1, a1, float_mean1, model_mean1)
    f_out2    = writefile(file_out2, args.var, count2, sigma2, corr_coeff2, bias_val2, b2, a2, float_mean2, model_mean2)

print('Calculation successfully computed.')
