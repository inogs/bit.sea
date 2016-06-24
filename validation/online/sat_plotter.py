from sat_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
import pylab as pl
import matplotlib.dates as mdates
from basins import V2 as OGS

TI_V1 = TimeInterval('20150407','20160607','%Y%m%d')
TI_V2 = TimeInterval('20160412','20160621','%Y%m%d')
ARCHIVE_DIR      ="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2-dev"
ARCHIVE_V1       ="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V4"
OUTFIG_DIR       ="/pico/home/userexternal/gcossari/COPERNICUS/CATENA/FIG_VALIDATION_ONLINE"

F0v1    = timelistcontainer(TI_V1,ARCHIVE_V1,'f0')
F1v1    = timelistcontainer(TI_V1,ARCHIVE_V1,'f1')
F2v1    = timelistcontainer(TI_V1,ARCHIVE_V1,'f2')

F0v2 = timelistcontainer(TI_V2,ARCHIVE_DIR      ,'f0')
F1v2 = timelistcontainer(TI_V2,ARCHIVE_DIR      ,'f1')
F2v2 = timelistcontainer(TI_V2,ARCHIVE_DIR      ,'f2')

coast='open_sea'
#for sub in [OGS.alb]:
for sub in OGS.P:
    print sub.name
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    # BIAS
    # f0  v1 
    times, f0v1 = F0v1.plotdata(F0v1.bias, sub.name, coast)
    ax1.plot(times,f0v1,'-b')
    times, f1v1 = F1v1.plotdata(F1v1.bias, sub.name, coast)
    ax1.plot(times,f1v1,'--b')
    times, f2v1 = F2v1.plotdata(F2v1.bias, sub.name, coast)
    ax1.plot(times,f2v1,':b')
    # f0  v2  
    times, f0= F0v2.plotdata(F0v2.bias, sub.name,coast)
    ax1.plot(times,f0,'-k',label='$1^{st}$ day of forecast (misfit DA d+7)')
    # f1  v2 
    times, f1= F1v2.plotdata(F1v2.bias, sub.name,coast)
    ax1.plot(times,f1,'--k',label='$2^{nd}$ day of forecast')
    # f2  v2  
    times, f2= F2v2.plotdata(F2v2.bias, sub.name,coast)
    ax1.plot(times,f2,':k',label='$3^{rd}$ day of forecast')
    ax1.set_title(sub.extended_name + " (" + coast + ")",fontsize=14)
    ax1.set_ylabel('BIAS $mg/m^3$', fontsize=14)
    ax1.legend(bbox_to_anchor=(.31,1.15), fontsize = 12)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax1.grid(True)    
    #RMS

    times, f0v1 = F0v1.plotdata(F0v1.rmse, sub.name, coast)
    ax2.plot(times,f0v1,'-b',label='V1')
    times, f1v1 = F1v1.plotdata(F1v1.rmse, sub.name, coast)
    ax2.plot(times,f1v1,'--b')
    times, f2v1 = F2v1.plotdata(F2v1.rmse, sub.name, coast)
    ax2.plot(times,f2v1,':b')
    # f0  v2  
    times, f0= F0v2.plotdata(F0v2.rmse, sub.name,coast)
    ax2.plot(times,f0,'-k',label='V2') 
    # f1  v2 
    times, f1= F1v2.plotdata(F1v2.rmse, sub.name,coast)
    ax2.plot(times,f1,'--k')
    # f2  v2  
    times, f2= F2v2.plotdata(F2v2.rmse, sub.name,coast)
    ax2.plot(times,f2,':k')
    ax2.set_ylabel('RMS $mg/m^3$',fontsize=14)
    ax2.legend(bbox_to_anchor=(.15,1.05), fontsize = 14)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax2.grid(True)

    # n. of valid points
    times, f0v1 = F0v1.plotdata(F0v1.number, sub.name, coast)
    ax3.plot(times,f0v1,'-b')
    times, f1v1 = F1v1.plotdata(F1v1.number, sub.name, coast)
    ax3.plot(times,f1v1,'--b')
    times, f2v1 = F2v1.plotdata(F2v1.number, sub.name, coast)
    ax3.plot(times,f2v1,':b')
    # f0  v2  
    times, f0= F0v2.plotdata(F0v2.number, sub.name,coast)
    ax3.plot(times,f0,'-k')
    # f1  v2 
    times, f1= F1v2.plotdata(F1v2.number, sub.name,coast)
    ax3.plot(times,f1,'--k')
    # f2  v2  
    times, f2= F2v2.plotdata(F2v2.number, sub.name,coast)
    ax3.plot(times,f2,':k')
    ax3.set_ylabel('# of points',fontsize=14)
    ax3.xaxis.set_major_locator(mdates.MonthLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))

    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    pl.show(block=False)
    nomefile=OUTFIG_DIR + "/NRTvalidation_chlsup_" + sub.name + "_" + coast + ".png"
    fig.savefig(nomefile)
    pl.close(fig)


##############















