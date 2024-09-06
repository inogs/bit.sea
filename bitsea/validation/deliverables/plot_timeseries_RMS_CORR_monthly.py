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

if (1 == 0): 
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
    pl.ylim(-0.5, 0.5)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
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


#from commons.season import season
#S=season()
#S.setseasons(["0101", "0501", "0601", "1001"], ["winter","spring","summer","fall"])
from commons import timerequestors
from commons.Timelist import TimeInterval, TimeList
TL=TimeList(TIMES)
from commons.utils import writetable

#iSeas=0 # JAN-APR
#CLIM_REQ=timerequestors.Clim_season(iSeas,S)
#ii,w=TL.select(CLIM_REQ)
#RMS__win = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
#BIAS_win = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
#RMSL_win = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
#BIASLwin = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

#MEAN_MOD_win = np.nanmean(MODEL_MEAN[ii,:],axis=0)
#MEAN_REF_win = np.nanmean(SAT___MEAN[ii,:],axis=0)
#
#STD_MOD_win = np.nanmean(MODEL_STD[ii,:],axis=0)
#STD_REF_win = np.nanmean(SAT___STD[ii,:],axis=0)
#CORR_win    = np.nanmean(BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

#iSeas=2 # JUN-SEP
#CLIM_REQ=timerequestors.Clim_season(iSeas,S)

RMS__m = np.zeros((12,nSUB),np.float32)
BIAS_m = np.zeros((12,nSUB),np.float32)
RMSL_m = np.zeros((12,nSUB),np.float32)
BIASLm = np.zeros((12,nSUB),np.float32)
MEAN_MOD_m = np.zeros((12,nSUB),np.float32)
MEAN_REF_m = np.zeros((12,nSUB),np.float32)
STD_MOD_m = np.zeros((12,nSUB),np.float32)
STD_REF_m = np.zeros((12,nSUB),np.float32)
CORR_m = np.zeros((12,nSUB),np.float32)

mat = np.zeros((12,nSUB,9),np.float32)

for imonth in range(12):
    req=timerequestors.Monthly_req(2019,imonth+1) 
    ii, w = TL.select(req)
    RMS__m[imonth,:] = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
    BIAS_m[imonth,:] = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
    RMSL_m[imonth,:] = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
    BIASLm[imonth,:]= np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)
    
    MEAN_MOD_m[imonth,:] = np.nanmean(MODEL_MEAN[ii,:],axis=0)
    MEAN_REF_m[imonth,:] = np.nanmean(SAT___MEAN[ii,:],axis=0)
    STD_MOD_m[imonth,:]  = np.nanmean(MODEL_STD[ii,:],axis=0)
    STD_REF_m[imonth,:]  = np.nanmean(SAT___STD[ii,:],axis=0)
    CORR_m[imonth,:]     = np.nanmean(BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

    mat[imonth,:,0] = RMS__m[imonth,:]
    mat[imonth,:,1] = BIAS_m[imonth,:]
    mat[imonth,:,2] = RMSL_m[imonth,:]
    mat[imonth,:,3] = BIASLm[imonth,:]
    mat[imonth,:,4] = STD_MOD_m[imonth,:]
    mat[imonth,:,5] = STD_REF_m[imonth,:]
    mat[imonth,:,6] = CORR_m[imonth,:]
#----
    mat[imonth,:,7] = MEAN_MOD_m[imonth,:]
    mat[imonth,:,8] = MEAN_REF_m[imonth,:]


for isub, sub in enumerate(OGS.P.basin_list):
    if (sub.name == "adr1"):
        outfiletable = args.outdir+"/"+"table4.1_month.dat"
#rows_names=[sub.name for sub in OGS.P.basin_list]
        rows_names=[imonth for imonth in np.arange(12)+1]
#column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum']
#column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum','STD_MODwin','STD_SATwin','STD_MODsum','STD_SATsum','CORRwin','CORRsum']
        column_names=['RMS','BIAS', 'RMSL','BIASL','STD_MOD','STD_REF','CORR','MEAN_MOD','MEAN_SAT']
        writetable(outfiletable, mat[:,isub,:], rows_names, column_names, fmt='%5.3f\t')

