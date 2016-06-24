from sat_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
import pylab as pl
import matplotlib.dates as mdates
from basins import V2 as OGS

TI = TimeInterval('20160401','20160601','%Y%m%d')
ARCHIVE_DIR_WRONG="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2"
ARCHIVE_DIR      ="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2-dev"
OUTFIG_DIR       ="/pico/home/userexternal/gcossari/COPERNICUS/CATENA/FIG_VALIDATION_ONLINE"

F0v2wrong = timelistcontainer(TI,ARCHIVE_DIR_WRONG,'f0')
F1v2wrong = timelistcontainer(TI,ARCHIVE_DIR_WRONG,'f1')
F2v2wrong = timelistcontainer(TI,ARCHIVE_DIR_WRONG,'f2')



F0v2 = timelistcontainer(TI,ARCHIVE_DIR      ,'f0')
F1v2 = timelistcontainer(TI,ARCHIVE_DIR      ,'f1')
F2v2 = timelistcontainer(TI,ARCHIVE_DIR      ,'f2')

for sub in [OGS.alb]:
    print sub.name
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    # BIAS
    # f0  v2 wrong 
    times, f0_w= F0v2wrong.plotdata(F0v2wrong.bias, sub.name,'open_sea')
    ax1.plot(times,f0_w,'-r',label='Wr: $1^{st}$ day of forecast (misfit DA d+7)')
    # f1  v2 wrong 
    times, f1_w= F1v2wrong.plotdata(F1v2wrong.bias, sub.name,'open_sea')
    ax1.plot(times,f1_w,'--r',label='Wr: $2^{nd}$ day of forecast')
    # f2  v2 wrong 
    times, f2_w= F2v2wrong.plotdata(F2v2wrong.bias, sub.name,'open_sea')
    ax1.plot(times,f2_w,':r',label='Wr: $3^{rd}$ day of forecast')
    # f0  v2  
    times, f0= F0v2.plotdata(F0v2.bias, sub.name,'open_sea')
    ax1.plot(times,f0,'-k',label='$1^{st}$ day of forecast (misfit DA d+7)')
    # f1  v2 
    times, f1= F1v2.plotdata(F1v2.bias, sub.name,'open_sea')
    ax1.plot(times,f1,'--k',label='$2^{nd}$ day of forecast')
    # f2  v2  
    times, f2= F2v2.plotdata(F2v2.bias, sub.name,'open_sea')
    ax1.plot(times,f2,':k',label='$3^{rd}$ day of forecast')
    ax1.set_title(sub.extended_name,fontsize=14)
    ax1.set_ylabel('BIAS $mg/m^3$', fontsize=14)
    ax1.legend(bbox_to_anchor=(1.05,1.15), fontsize = 12)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
        
    #RMS
    # f0  v2 wrong 
    times, f0_w= F0v2wrong.plotdata(F0v2wrong.rmse, sub.name,'open_sea')
    ax2.plot(times,f0_w,'-r')
    # f1  v2 wrong 
    times, f1_w= F1v2wrong.plotdata(F1v2wrong.rmse, sub.name,'open_sea')
    ax2.plot(times,f1_w,'--r') 
    # f2  v2 wrong 
    times, f2_w= F2v2wrong.plotdata(F2v2wrong.rmse, sub.name,'open_sea')
    ax2.plot(times,f2_w,':r')
    # f0  v2  
    times, f0= F0v2.plotdata(F0v2.rmse, sub.name,'open_sea')
    ax2.plot(times,f0,'-k') 
    # f1  v2 
    times, f1= F1v2.plotdata(F1v2.rmse, sub.name,'open_sea')
    ax2.plot(times,f1,'--k')
    # f2  v2  
    times, f2= F2v2.plotdata(F2v2.rmse, sub.name,'open_sea')
    ax2.plot(times,f2,':k')
    ax2.set_ylabel('RMS $mg/m^3$',fontsize=14)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))

    # n. of valid points
    # f0  v2 wrong 
    times, f0_w= F0v2wrong.plotdata(F0v2wrong.number, sub.name,'open_sea')
    ax3.plot(times,f0_w,'-r')
    # f1  v2 wrong 
    times, f1_w= F1v2wrong.plotdata(F1v2wrong.number, sub.name,'open_sea')
    ax3.plot(times,f1_w,'--r')
    # f2  v2 wrong 
    times, f2_w= F2v2wrong.plotdata(F2v2wrong.number, sub.name,'open_sea')
    ax3.plot(times,f2_w,':r')
    # f0  v2  
    times, f0= F0v2.plotdata(F0v2.number, sub.name,'open_sea')
    ax3.plot(times,f0,'-k')
    # f1  v2 
    times, f1= F1v2.plotdata(F1v2.number, sub.name,'open_sea')
    ax3.plot(times,f1,'--k')
    # f2  v2  
    times, f2= F2v2.plotdata(F2v2.number, sub.name,'open_sea')
    ax3.plot(times,f2,':k')
    ax3.set_ylabel('# of points',fontsize=14)
    ax3.xaxis.set_major_locator(mdates.MonthLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))

    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    pl.show(block=False)
    nomefile=OUTFIG_DIR + "/Fig_NRTvalidation_sat_" + sub.name + ".png"
    fig.savefig(nomefile)
    pl.close(fig)


##############















