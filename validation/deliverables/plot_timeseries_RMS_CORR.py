# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.3 and table IV.1


import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import sys
import numpy as np
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )

    parser.add_argument(   '--inputfile', '-i',
                            type = str,
                            required = True,
                            default = 'export_data_ScMYValidation_plan.pkl',
                            help = 'Input pickle file')


    return parser.parse_args()

args = argument()



fid = open(args.inputfile)
LIST = pickle.load(fid)
fid.close()

TIMES                          = LIST[0]
BGC_CLASS4_CHL_RMS_SURF_BASIN  = LIST[1]
BGC_CLASS4_CHL_BIAS_SURF_BASIN = LIST[2]
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = LIST[5]
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= LIST[6]

MODEL_MEAN=LIST[3]
SAT___MEAN=LIST[4]
MODEL_STD = LIST[7]
SAT___STD = LIST[8]
BGC_CLASS4_CHL_CORR_SURF_BASIN= LIST[9]

from basins import V2 as OGS

nSUB = len(OGS.P.basin_list)

for isub,sub in enumerate(OGS.P):
#  if (isub != 17):
    print sub.name
    fig, ax = pl.subplots()
    ax.plot(TIMES,BGC_CLASS4_CHL_RMS_SURF_BASIN[:,isub],'-k',label='RMS')
    ax.plot(TIMES,BGC_CLASS4_CHL_BIAS_SURF_BASIN[:,isub],'-b',label='Bias')
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(-0.3, 0.3)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=9)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=10)
#    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename=args.outdir+"/"+'chl-RMS-BIAS_' + sub.name + ".png"
    pl.savefig(outfilename)
    pl.close(fig)


from commons.season import season
S=season()
S.setseasons(["0101", "0501", "0601", "1001"], ["winter","spring","summer","fall"])
from commons import timerequestors
from commons.Timelist import TimeInterval, TimeList
TL=TimeList(TIMES)
from commons.utils import writetable

iSeas=0 # JAN-APR
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__win = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_win = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_win = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLwin = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

MEAN_MOD_win = np.nanmean(MODEL_MEAN[ii,:],axis=0)
MEAN_REF_win = np.nanmean(SAT___MEAN[ii,:],axis=0)

STD_MOD_win = np.nanmean(MODEL_STD[ii,:],axis=0)
STD_REF_win = np.nanmean(SAT___STD[ii,:],axis=0)
CORR_win    = np.nanmean(BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

iSeas=2 # JUN-SEP
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__sum = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_sum = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_sum = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLsum = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

MEAN_MOD_sum = np.nanmean(MODEL_MEAN[ii,:],axis=0)
MEAN_REF_sum = np.nanmean(SAT___MEAN[ii,:],axis=0)

STD_MOD_sum = np.nanmean(MODEL_STD[ii,:],axis=0)
STD_REF_sum = np.nanmean(SAT___STD[ii,:],axis=0)
CORR_sum    = np.nanmean(BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

#mat = np.zeros((nSUB,8),np.float32)
#mat = np.zeros((nSUB,14),np.float32)
mat = np.zeros((nSUB,18),np.float32)

mat[:,0] = RMS__win
mat[:,1] = RMS__sum
mat[:,2] = BIAS_win
mat[:,3] = BIAS_sum
mat[:,4] = RMSL_win
mat[:,5] = RMSL_sum
mat[:,6] = BIASLwin
mat[:,7] = BIASLsum
mat[:,8] = STD_MOD_win
mat[:,9] = STD_REF_win
mat[:,10] = STD_MOD_sum
mat[:,11] = STD_REF_sum
mat[:,12] = CORR_win
mat[:,13] = CORR_sum
#----
mat[:,14] = MEAN_MOD_win
mat[:,15] = MEAN_REF_win
mat[:,16] = MEAN_MOD_sum
mat[:,17] = MEAN_REF_sum

outfiletable = args.outdir+"/"+"table4.1.dat"
rows_names=[sub.name for sub in OGS.P.basin_list]
#column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum']
#column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum','STD_MODwin','STD_SATwin','STD_MODsum','STD_SATsum','CORRwin','CORRsum']
column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum','STD_MODwin','STD_SATwin','STD_MODsum','STD_SATsum','CORRwin','CORRsum','MEAN_MODwin','MEAN_SATwin','MEAN_MODsum','MEAN_SATsum']
writetable(outfiletable, mat, rows_names, column_names, fmt='%5.3f\t')

