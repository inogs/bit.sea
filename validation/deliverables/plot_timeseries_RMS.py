# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.3


import pickle
import pylab as pl
import matplotlib.dates as mdates
import sys
import numpy as np

fid = open('export_data_ScMYValidation_plan.pkl')
LIST = pickle.load(fid)
fid.close()
TIMES                          = LIST[0]
BGC_CLASS4_CHL_RMS_SURF_BASIN  = LIST[1]
BGC_CLASS4_CHL_BIAS_SURF_BASIN = LIST[2]
MODEL_MEAN                     = LIST[3]
SAT___MEAN                     = LIST[4]
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = LIST[5]
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= LIST[6]

#surf_layer = Layer(0,10)

from basins import OGS 
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
    pl.setp(xlabels, rotation=30)
#    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    #outfilename='chl-RMS-BIAS_' + sub.name + ".png"
    #pl.savefig(outfilename)
    #sys.exit()

start_year=1999
end___year=2012
nyr=end___year-start_year+1
RMS__sum=np.zeros((nyr,11),np.float32)
RMS__win=np.zeros((nyr,11),np.float32)
BIAS_sum=np.zeros((nyr,11),np.float32)
BIAS_win=np.zeros((nyr,11),np.float32)
RMSL_sum=np.zeros((nyr,11),np.float32)
RMSL_win=np.zeros((nyr,11),np.float32)
BIASLsum=np.zeros((nyr,11),np.float32)
BIASLwin=np.zeros((nyr,11),np.float32)
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
    RMS__win[n,:]=LIST[1][w1:w2,:].mean(axis=0)
    RMS__sum[n,:]=LIST[1][s1:s2,:].mean(axis=0)
    BIAS_win[n,:]=LIST[2][w1:w2,:].mean(axis=0)
    BIAS_sum[n,:]=LIST[2][s1:s2,:].mean(axis=0)
    RMSL_win[n,:]=LIST[5][w1:w2,:].mean(axis=0)
    RMSL_sum[n,:]=LIST[5][s1:s2,:].mean(axis=0)
    BIASLwin[n,:]=LIST[6][w1:w2,:].mean(axis=0)
    BIASLsum[n,:]=LIST[6][s1:s2,:].mean(axis=0)

#TABELLA QUID IV.1    
print 'RMS_win : ',RMS__win[:,:].mean(axis=0)
print 'RMS_sum : ',RMS__sum[:,:].mean(axis=0)
print 'BIAS_win: ',BIAS_win[:,:].mean(axis=0)
print 'BIAS_sum: ',BIAS_sum[:,:].mean(axis=0)
print 'RMSLwin : ',RMSL_win[:,:].mean(axis=0)
print 'RMSLsum : ',RMSL_sum[:,:].mean(axis=0)
print 'BIASLwin: ',BIASLwin[:,:].mean(axis=0)
print 'BIASLsum: ',BIASLsum[:,:].mean(axis=0)







