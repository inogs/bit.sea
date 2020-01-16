import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates time series png files for mdeeaf web site.
    Each file has 3 subplots: bias, rms and number of points.
    The time window is 18 months.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--date','-d',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')

    parser.add_argument(   '--archivedir','-a',
                                type = str,
                                required = True,
                                help = 'chain archive directory')
    parser.add_argument(   '--validation_dir','-v',
                                type = str,
                                required = True,
                                help = 'wrkdir/2/POSTPROC/AVE_FREQ_1/online_validation')
#    parser.add_argument(   '--previous_archive','-p',
#                                type = str,
#                                required = True,
#                                help = 'previous chain archive directory, taken from static-data')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import FormatStrFormatter

from sat_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
from commons.utils import addsep

import matplotlib.pyplot as pl
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from datetime import datetime
from basins import V2 as OGS
from commons.utils import writetable

nSub=len(OGS.P.basin_list)
Graphic_DeltaT = relativedelta(months=18)
datestart = datetime.strptime(args.date,'%Y%m%d') -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")
##TI_V1 = TimeInterval(timestart,'20160607','%Y%m%d')
##TI_V2 = TimeInterval('20160412', args.date,'%Y%m%d')
#TI_V2 = TimeInterval(timestart,'20171107','%Y%m%d')
TI_V3 = TimeInterval('20180501', args.date,'%Y%m%d')
ARCHIVE_DIR      = addsep(args.archivedir) #"/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2-dev"
#ARCHIVE_PREV     = addsep(args.previous_archive)   #r"/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V4"
VALID_WRKDIR     = addsep(args.validation_dir)
OUTFIG_DIR       = addsep(args.outdir) # "/pico/home/userexternal/gcossari/COPERNICUS/CATENA/FIG_VALIDATION_ONLINE"
var = "chlsup"
##F0v1    = timelistcontainer(TI_V1,ARCHIVE_PREV,'f0', postfix_dir="")
##F1v1    = timelistcontainer(TI_V1,ARCHIVE_PREV,'f1', postfix_dir="")
##F2v1    = timelistcontainer(TI_V1,ARCHIVE_PREV,'f2', postfix_dir="")
#print TI_V2
print TI_V3
#print ARCHIVE_PREV
#F0v2    = timelistcontainer(TI_V2,ARCHIVE_PREV,'f0', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")
#F1v2    = timelistcontainer(TI_V2,ARCHIVE_PREV,'f1', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")
#F2v2    = timelistcontainer(TI_V2,ARCHIVE_PREV,'f2', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")

#F0v3 = timelistcontainer(TI_V3,ARCHIVE_DIR, 'f0', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")
F0v3 = timelistcontainer(TI_V3,"/gpfs/scratch/userexternal/lfeudale/work/OPA/SAT_f0_DAILY/tmp/" , 'f0' , postfix_dir="")
F1v3 = timelistcontainer(TI_V3,ARCHIVE_DIR, 'f1', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")
F2v3 = timelistcontainer(TI_V3,ARCHIVE_DIR, 'f2', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")

F0v3.append_dir(VALID_WRKDIR)
F1v3.append_dir(VALID_WRKDIR)
F1v3.append_dir(VALID_WRKDIR)

coast='open_sea'
column_names=["f0 MEAN","f1 MEAN","f2 MEAN"]
row_names   =[sub.name for sub in OGS.P.basin_list]
RMSD_mean = np.zeros((nSub,3),np.float32)
#STAT_mean = np.zeros((nSub,21),np.float32)
STAT_mean = np.zeros((nSub,18),np.float32)

#EAN_BIAS_w = [ ]
EAN_BIAS_s = [0.02, -0.01, -0.01, -0.01, -0.01, 0.005, -0.02, -0.01, -0.01, 0.005, 0.005, -0.01, 0.005, 0.005, 0.005, 0.005, 0.005]
EAN_RMSD_w = [0.18, 0.07, 0.05, 0.10, 0.05, 0.04, 0.03, 0.04, 0.03, 0.02, 0.01, 0.03, 0.02, 0.02, 0.01, 0.02, 0.06 ]
EAN_RMSD_s = [0.09, 0.01, 0.01, 0.02, 0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.005, 0.01, 0.005, 0.005, 0.005, 0.01, 0.02]
#for sub in [OGS.alb]:
for isub, sub in enumerate(OGS.P):
    print sub.name
    outfile=OUTFIG_DIR + "NRTvalidation_chlsup_" + sub.name + "_" + coast + "_RMSD_f0daily.png"
#    #matplotlib.rc('xtick', labelsize=12)
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    fig, ax2 = pl.subplots(1, sharex=True, figsize=(15,10))
    fig.text(0.04,0.93,sub.name,fontsize=34,color="r")
    # BIAS
    # v3  
#    times, f0= F0v3.plotdata(F0v3.bias, sub.name,coast) ; ax1.plot(times,f0,'-k' ,label='1$^{st}$ day of forecast (misfit DA d+7)')
#    ax1.plot(times,EAN_BIAS[isub]*np.ones(len(times)),'-.r')
#    times, f1= F1v3.plotdata(F1v3.bias, sub.name,coast);  ax1.plot(times,f1,'--k',label='2$^{nd}$ day of forecast')  
#    times, f2= F2v3.plotdata(F2v3.bias, sub.name,coast);  ax1.plot(times,f2,':k' ,label='3$^{rd}$ day of forecast')
#    
#    ax1.set_title(sub.extended_name ,fontsize=14)
#    ax1.set_ylabel('BIAS mg/m$^3$', fontsize=14)
#    max_val1=max(list(np.abs(f0[~np.isnan(f0)])) + list(np.abs(f1[~np.isnan(f1)])) + list(np.abs(f2[~np.isnan(f2)])) )
#    ax1.set_ylim(-1.2*max_val1,1.2*max_val1)
#    ax1.tick_params(axis='both', labelsize=14)
#    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#    ax1.grid(True)
    
    
    #RMS
    # fv3  
    times, f0= F0v3.plotdata(F0v3.rmse, sub.name,coast) ; ax2.plot(times,f0,'-k',label='V3') 
    print "times0"
    print times , f0
    EAN_RMSD=np.ndarray((len(times),nSub), dtype = np.float32)
    for it, t in enumerate(times):
        print times[it].month
        if (( t.month <= 5 ) | ( t.month >= 10 )): 
           EAN_RMSD[it,isub] = EAN_RMSD_w[isub]
        else:
           EAN_RMSD[it,isub] = EAN_RMSD_s[isub]
    ax2.plot(times,EAN_RMSD[:,isub]*np.ones(len(times)),'-.r')
    print "again"
    print times
    times, f1= F1v3.plotdata(F1v3.rmse, sub.name,coast) ; ax2.plot(times,f1,'--k') 
    print "times1"
    print times
    times, f2= F2v3.plotdata(F2v3.rmse, sub.name,coast) ; ax2.plot(times,f2,':k')
    print "times2"
    print times
    max_val2=max(list(np.abs(f0[~np.isnan(f0)])) + list(np.abs(f1[~np.isnan(f1)])) + list(np.abs(f2[~np.isnan(f2)])) )
    ax2.set_ylim(0,1.2*max_val2)  
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax2.tick_params(axis='both', labelsize=28)
    ax2.set_ylabel('RMSD mg/m$^3$',fontsize=28)
    ax2.grid(True)

    RMSD_mean[:,0]=np.nanmean(F0v3.rmse[:,:,1],axis=0)
    RMSD_mean[:,1]=np.nanmean(F1v3.rmse[:,:,1],axis=0)
    RMSD_mean[:,2]=np.nanmean(F2v3.rmse[:,:,1],axis=0)

    STAT_mean[:,0]=np.nanmean(F0v3.rmse[:,:,1],axis=0)
    STAT_mean[:,1]=np.nanmean(F1v3.rmse[:,:,1],axis=0)
    STAT_mean[:,2]=np.nanmean(F2v3.rmse[:,:,1],axis=0)
                   
    STAT_mean[:,3]=np.nanmean(F0v3.model[:,:,1],axis=0)
    STAT_mean[:,4]=np.nanmean(F1v3.model[:,:,1],axis=0)
    STAT_mean[:,5]=np.nanmean(F2v3.model[:,:,1],axis=0)
    STAT_mean[:,6]=np.nanmean(F0v3.sat[:,:,1],axis=0)
    STAT_mean[:,7]=np.nanmean(F1v3.sat[:,:,1],axis=0)
    STAT_mean[:,8]=np.nanmean(F2v3.sat[:,:,1],axis=0)
    STAT_mean[:,9]=np.nanstd(F0v3.model[:,:,1],axis=0)
    STAT_mean[:,10]=np.nanstd(F1v3.model[:,:,1],axis=0)
    STAT_mean[:,11]=np.nanstd(F2v3.model[:,:,1],axis=0)
    STAT_mean[:,12]=np.nanstd(F0v3.sat[:,:,1],axis=0)
    STAT_mean[:,13]=np.nanstd(F1v3.sat[:,:,1],axis=0)
    STAT_mean[:,14]=np.nanstd(F2v3.sat[:,:,1],axis=0)

    STAT_mean[:,15]=np.nanmean(F0v3.bias[:,:,1],axis=0)
    STAT_mean[:,16]=np.nanmean(F1v3.bias[:,:,1],axis=0)
    STAT_mean[:,17]=np.nanmean(F2v3.bias[:,:,1],axis=0)
#    STAT_mean[:,18]=np.corrcoef((F0v3.sat[:,:,1],F0v3.model[:,:,1]),axis=0)
#    STAT_mean[:,19]=np.corrcoef((F1v3.sat[:,:,1],F1v3.model[:,:,1]),axis=0)
#    STAT_mean[:,20]=np.corrcoef((F2v3.sat[:,:,1],F2v3.model[:,:,1]),axis=0)
    
    # v3  
 #   times, f0= F0v3.plotdata(F0v3.number, sub.name,coast)     ; ax3.plot(times,f0,'-k') 
 #   ax3.xaxis_date()
 #   ax3.set_yticklabels(ax3.get_yticks() ,fontsize=14)
 #   ax3.grid(True)
 #   ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))
 #   pl.tick_params(labelsize = 14)
 #   times, f1= F1v3.plotdata(F1v3.number, sub.name,coast)     ; ax3.plot(times,f1,'--k')
 #   times, f2= F2v3.plotdata(F2v3.number, sub.name,coast)     ; ax3.plot(times,f2,':k')
 #   xlabels = ax3.get_xticklabels()
 #   pl.setp(xlabels, rotation=30) 
 #   ax3.set_ylabel('# of points',fontsize=14)

    column_names_STAT=["f0_RMDS_MEAN","f1_RMSD_MEAN","f2_RMSD_MEAN","f0_MOD_MEAN","f1_MOD_MEAN","f2_MOD_MEAN","f0_SAT_MEAN","f1_SAT_MEAN","f2_SAT_MEAN","f0_MOD_STD","f1_MOD_STD","f2_MOD_STD","f0_SAT_STD","f1_SAT_STD","f2_SAT_STD","f0_BIAS","f1_BIAS","f2_BIAS"]
    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.savefig(outfile)
    pl.close(fig)
    writetable(OUTFIG_DIR + 'NRT_fc_' + var + '_RMSD_f0daily.txt',RMSD_mean,row_names, column_names)
    writetable(OUTFIG_DIR + 'NRT_fc_' + var + '_STAT_f0daily.txt',STAT_mean,row_names, column_names_STAT)
