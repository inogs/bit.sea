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

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')

from sat_timeseries import timelistcontainer
from commons.time_interval import TimeInterval
from commons.utils import addsep

import pylab as pl
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from datetime import datetime
from basins import V2 as OGS

Graphic_DeltaT = relativedelta(months=18)
datestart = datetime.strptime(args.date,'%Y%m%d') -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")


TI_V3 = TimeInterval('20190405', args.date,'%Y%m%d')
ARCHIVE_DIR      = addsep(args.archivedir) #"/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/V2-dev"
VALID_WRKDIR     = addsep(args.validation_dir)
OUTFIG_DIR       = addsep(args.outdir) # "/pico/home/userexternal/gcossari/COPERNICUS/CATENA/FIG_VALIDATION_ONLINE"

print TI_V3

postfix="POSTPROC/AVE_FREQ_1/validation/Sat/"
postfix=""
F0v3 = timelistcontainer(TI_V3,ARCHIVE_DIR, 'f0', postfix_dir=postfix)
F1v3 = timelistcontainer(TI_V3,ARCHIVE_DIR, 'f1', postfix_dir=postfix)
F2v3 = timelistcontainer(TI_V3,ARCHIVE_DIR, 'f2', postfix_dir=postfix)

F0v3.append_dir(VALID_WRKDIR)
F1v3.append_dir(VALID_WRKDIR)
F1v3.append_dir(VALID_WRKDIR)

coast='open_sea'
#for sub in [OGS.alb]:
for sub in OGS.P:
    print sub.name
    outfile=OUTFIG_DIR + "NRTvalidation_chlsup_" + sub.name + "_" + coast + ".png"
    #matplotlib.rc('xtick', labelsize=12)
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    # BIAS
    # v3  
    times, f0= F0v3.plotdata(F0v3.bias, sub.name,coast) ; ax1.plot(times,f0,'-k' ,label='1$^{st}$ day of forecast (misfit DA d+7)')
    times, f1= F1v3.plotdata(F1v3.bias, sub.name,coast);  ax1.plot(times,f1,'--k',label='2$^{nd}$ day of forecast')  
    times, f2= F2v3.plotdata(F2v3.bias, sub.name,coast);  ax1.plot(times,f2,':k' ,label='3$^{rd}$ day of forecast')
    
    ax1.set_title(sub.extended_name + " (" + coast + ")",fontsize=14)
    ax1.set_ylabel('BIAS mg/m$^3$', fontsize=14)
#HERE
#    ax1.legend(bbox_to_anchor=(.31,1.15), fontsize = 12)
    ax1.legend(bbox_to_anchor=(0.37,1), bbox_transform=pl.gcf().transFigure, fontsize = 12)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.grid(True)
    
    #RMS
    # fv3  
    times, f0= F0v3.plotdata(F0v3.rmse, sub.name,coast) ; ax2.plot(times,f0,'-k',label='V3') 
    times, f1= F1v3.plotdata(F1v3.rmse, sub.name,coast) ; ax2.plot(times,f1,'--k') 
    times, f2= F2v3.plotdata(F2v3.rmse, sub.name,coast) ; ax2.plot(times,f2,':k')
    
    ax2.set_ylabel('RMS mg/m$^3$',fontsize=14)
    ax2.legend(bbox_to_anchor=(.15,1.05), fontsize = 14)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.grid(True)


    # n. of valid points
    # v3  
    times, f0= F0v3.plotdata(F0v3.number, sub.name,coast)     ; ax3.plot(times,f0,'-k') 
    times, f1= F1v3.plotdata(F1v3.number, sub.name,coast)     ; ax3.plot(times,f1,'--k')
    times, f2= F2v3.plotdata(F2v3.number, sub.name,coast)     ; ax3.plot(times,f2,':k')
    

    ax3.set_ylabel('# of points',fontsize=14)
    ax3.xaxis.set_major_locator(mdates.MonthLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%b-%y"))

    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.savefig(outfile)
    pl.close(fig)
    import sys
#    sys.exit()
