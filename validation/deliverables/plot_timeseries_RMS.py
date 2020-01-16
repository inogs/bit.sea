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


from basins import V2 as OGS

nSUB = len(OGS.P.basin_list)

for isub,sub in enumerate(OGS.P):
    print sub.name
    fig, ax = pl.subplots()
    ax.plot(TIMES,BGC_CLASS4_CHL_RMS_SURF_BASIN[:,isub],'-k',label='RMS')
    ax.hold(True)
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

iSeas=2 # JUN-SEP
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__sum = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_sum = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_sum = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLsum = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)
mat = np.zeros((nSUB,8),np.float32)

mat[:,0] = RMS__win
mat[:,1] = RMS__sum
mat[:,2] = BIAS_win
mat[:,3] = BIAS_sum
mat[:,4] = RMSL_win
mat[:,5] = RMSL_sum
mat[:,6] = BIASLwin
mat[:,7] = BIASLsum
outfiletable = args.outdir+"/"+"table4.1.dat"
rows_names=[sub.name for sub in OGS.P.basin_list]
column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum']
writetable(outfiletable, mat, rows_names, column_names, fmt='%5.3f\t')

import sys
sys.exit()


start_year=2015
end___year=2017
nyr=end___year-start_year+1
RMS__sum=np.zeros((nyr,nSUB),np.float32)
RMS__win=np.zeros((nyr,nSUB),np.float32)
BIAS_sum=np.zeros((nyr,nSUB),np.float32)
BIAS_win=np.zeros((nyr,nSUB),np.float32)
RMSL_sum=np.zeros((nyr,nSUB),np.float32)
RMSL_win=np.zeros((nyr,nSUB),np.float32)
BIASLsum=np.zeros((nyr,nSUB),np.float32)
BIASLwin=np.zeros((nyr,nSUB),np.float32)
#
# winter = JFMA
# summer = JJAS
#
for i in range(start_year,end___year+1):
    print i
    n=i-start_year
    w1=n*12
    w2=w1+4
    s1=w1+5
    s2=s1+4
    print w1,w2,s1,s2
    RMS__win[n,:]=BGC_CLASS4_CHL_RMS_SURF_BASIN[ w1:w2,:].mean(axis=0)
    RMS__sum[n,:]=BGC_CLASS4_CHL_RMS_SURF_BASIN[ s1:s2,:].mean(axis=0)
    BIAS_win[n,:]=BGC_CLASS4_CHL_BIAS_SURF_BASIN[w1:w2,:].mean(axis=0)
    BIAS_sum[n,:]=BGC_CLASS4_CHL_BIAS_SURF_BASIN[s1:s2,:].mean(axis=0)
    RMSL_win[n,:]=BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ w1:w2,:].mean(axis=0)
    RMSL_sum[n,:]=BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ s1:s2,:].mean(axis=0)
    BIASLwin[n,:]=BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[w1:w2,:].mean(axis=0)
    BIASLsum[n,:]=BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[s1:s2,:].mean(axis=0)

mat = np.zeros((8,nSUB))

mat[0,:] = RMS__win[:,:].mean(axis=0)
mat[1,:] = RMS__sum[:,:].mean(axis=0)
mat[2,:] = BIAS_win[:,:].mean(axis=0)
mat[3,:] = BIAS_sum[:,:].mean(axis=0)
mat[4,:] = RMSL_win[:,:].mean(axis=0)
mat[5,:] = RMSL_sum[:,:].mean(axis=0)
mat[6,:] = BIASLwin[:,:].mean(axis=0)
mat[7,:] = BIASLsum[:,:].mean(axis=0)

mat = mat.T
print '-----------------------------------------'
#TABELLA QUID IV.1
lines=[]
myformat = '%s ' + '%5.3g '*8 +"\n";
for isub,sub in enumerate(OGS.P):
    line = myformat %( sub.name, mat[isub,0], mat[isub,1], mat[isub,2],
                       mat[isub,3], mat[isub,4], mat[isub,5],
                       mat[isub,6],mat[isub,7])

    lines.append(line)

outfiletable = args.outdir+"/"+"table4.1.dat"
file = open(outfiletable,"w")
file.writelines(lines)
file.close()
# RMS__win[:,:].mean(axis=0)
# print 'RMS_sum : ',RMS__sum[:,:].mean(axis=0)
# print 'BIAS_win: ',BIAS_win[:,:].mean(axis=0)
# print 'BIAS_sum: ',BIAS_sum[:,:].mean(axis=0)
# print 'RMSLwin : ',RMSL_win[:,:].mean(axis=0)
# print 'RMSLsum : ',RMSL_sum[:,:].mean(axis=0)
# print 'BIASLwin: ',BIASLwin[:,:].mean(axis=0)
# print 'BIASLsum: ',BIASLsum[:,:].mean(axis=0)
