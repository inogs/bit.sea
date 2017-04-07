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
    parser.add_argument(   '--previous_archive','-p',
                                type = str,
                                required = True,
                                help = 'previous chain archive directory, taken from static-data')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

from sat_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
from commons.utils import addsep
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from datetime import datetime
from basins import V2 as OGS

Graphic_DeltaT = relativedelta(months=18)
datestart = datetime.strptime(args.date,'%Y%m%d') -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")
TI_V1 = TimeInterval(timestart,'20160607','%Y%m%d')
TI_V2 = TimeInterval('20160412', args.date,'%Y%m%d')
ARCHIVE_DIR      = addsep(args.archivedir) #"/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2-dev"
ARCHIVE_PREV     = addsep(args.previous_archive)   #r"/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V4"
VALID_WRKDIR     = addsep(args.validation_dir)
OUTFIG_DIR       = addsep(args.outdir) # "/pico/home/userexternal/gcossari/COPERNICUS/CATENA/FIG_VALIDATION_ONLINE"

F0v1    = timelistcontainer(TI_V1,ARCHIVE_PREV,'f0', postfix_dir="")
F1v1    = timelistcontainer(TI_V1,ARCHIVE_PREV,'f1', postfix_dir="")
F2v1    = timelistcontainer(TI_V1,ARCHIVE_PREV,'f2', postfix_dir="")

F0v2 = timelistcontainer(TI_V2,ARCHIVE_DIR, 'f0', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")
F1v2 = timelistcontainer(TI_V2,ARCHIVE_DIR, 'f1', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")
F2v2 = timelistcontainer(TI_V2,ARCHIVE_DIR, 'f2', postfix_dir="POSTPROC/AVE_FREQ_1/validation/Sat/")

F0v2.append_dir(VALID_WRKDIR)
F1v2.append_dir(VALID_WRKDIR)
F1v2.append_dir(VALID_WRKDIR)

coast='open_sea'
#for sub in [OGS.alb]:
for sub in OGS.P:
    print sub.name
    outfile=OUTFIG_DIR + "NRTvalidation_chlsup_" + sub.name + "_" + coast + ".png"
    #matplotlib.rc('xtick', labelsize=12)
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    # BIAS
    # v1 
    times, f0v1 = F0v1.plotdata(F0v1.bias, sub.name, coast); ax1.plot(times,f0v1,'-b')
    times, f1v1 = F1v1.plotdata(F1v1.bias, sub.name, coast); ax1.plot(times,f1v1,'--b')
    times, f2v1 = F2v1.plotdata(F2v1.bias, sub.name, coast); ax1.plot(times,f2v1,':b')
    # v2  
    times, f0= F0v2.plotdata(F0v2.bias, sub.name,coast) ; ax1.plot(times,f0,'-k' ,label='1$^{st}$ day of forecast (misfit DA d+7)')
    times, f1= F1v2.plotdata(F1v2.bias, sub.name,coast);  ax1.plot(times,f1,'--k',label='2$^{nd}$ day of forecast')  
    times, f2= F2v2.plotdata(F2v2.bias, sub.name,coast);  ax1.plot(times,f2,':k' ,label='3$^{rd}$ day of forecast')
    
    ax1.set_title(sub.extended_name + " (" + coast + ")",fontsize=14)
    ax1.set_ylabel('BIAS mg/m$^3$', fontsize=14)
    ax1.legend(bbox_to_anchor=(.31,1.15), fontsize = 12)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.grid(True)
    
    #RMS
    times, f0v1 = F0v1.plotdata(F0v1.rmse, sub.name, coast) ; ax2.plot(times,f0v1,'-b',label='V1')
    times, f1v1 = F1v1.plotdata(F1v1.rmse, sub.name, coast) ; ax2.plot(times,f1v1,'--b')
    times, f2v1 = F2v1.plotdata(F2v1.rmse, sub.name, coast) ; ax2.plot(times,f2v1,':b')
    # fv2  
    times, f0= F0v2.plotdata(F0v2.rmse, sub.name,coast) ; ax2.plot(times,f0,'-k',label='V2') 
    times, f1= F1v2.plotdata(F1v2.rmse, sub.name,coast) ; ax2.plot(times,f1,'--k') 
    times, f2= F2v2.plotdata(F2v2.rmse, sub.name,coast) ; ax2.plot(times,f2,':k')
    
    ax2.set_ylabel('RMS mg/m$^3$',fontsize=14)
    ax2.legend(bbox_to_anchor=(.15,1.05), fontsize = 14)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.grid(True)


    # n. of valid points
    times, f0v1 = F0v1.plotdata(F0v1.number, sub.name, coast) ; ax3.plot(times,f0v1,'-b')
    times, f1v1 = F1v1.plotdata(F1v1.number, sub.name, coast) ; ax3.plot(times,f1v1,'--b')
    times, f2v1 = F2v1.plotdata(F2v1.number, sub.name, coast) ; ax3.plot(times,f2v1,':b')
    # v2  
    times, f0= F0v2.plotdata(F0v2.number, sub.name,coast)     ; ax3.plot(times,f0,'-k') 
    times, f1= F1v2.plotdata(F1v2.number, sub.name,coast)     ; ax3.plot(times,f1,'--k')
    times, f2= F2v2.plotdata(F2v2.number, sub.name,coast)     ; ax3.plot(times,f2,':k')
    

    ax3.set_ylabel('# of points',fontsize=14)
    ax3.xaxis.set_major_locator(mdates.MonthLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%b-%y"))

    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.savefig(outfile)
    pl.close(fig)

