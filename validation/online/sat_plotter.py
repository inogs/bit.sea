import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates time series png files for mdeeaf web site.
    Each file has 3 subplots: bias, rms and number of points.
    The time window is 18 months.
    See https://matplotlib.org/3.1.0/gallery/color/named_colors.html for colors
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
import numpy as np
import pylab as pl
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from datetime import datetime
from basins import V2 as OGS

def avoid_diagonals_in_jumps(x,y):
    n = len(x)
    newY=[]
    newX=[]
    for i in range(n-1):
        if y[i]==y[i+1]:
           newY.append(y[i])
           newX.append(x[i])
        else:
            newX.append(x[i])
            newX.append(x[i])
            newY.append(y[i])
            newY.append(np.nan)
    newX.append(x[-1])
    newY.append(y[-1])
    return newX, newY


Graphic_DeltaT = relativedelta(months=18)
datestart = datetime.strptime(args.date,'%Y%m%d') -Graphic_DeltaT
timestart = datestart.strftime("%Y%m%d")


TI = TimeInterval('20190405', args.date,'%Y%m%d')
ARCHIVE_DIR      = addsep(args.archivedir)
VALID_WRKDIR     = addsep(args.validation_dir)
OUTFIG_DIR       = addsep(args.outdir)



prefix="Validation_f0_20190723_on_daily_Sat."
F0 = timelistcontainer(TI,ARCHIVE_DIR, 'Validation_f0_*on_daily_Sat*', prefix=prefix)
F1 = timelistcontainer(TI,ARCHIVE_DIR, 'Validation_f1_*on_daily_Sat*', prefix=prefix)
F2 = timelistcontainer(TI,ARCHIVE_DIR, 'Validation_f2_*on_daily_Sat*', prefix=prefix)

F0.append_dir(VALID_WRKDIR)
F1.append_dir(VALID_WRKDIR)
F1.append_dir(VALID_WRKDIR)

EAN_BIAS_s = [0.02, -0.01, -0.01, -0.01, -0.01, 0.005, -0.02, -0.01, -0.01, 0.005, 0.005, -0.01, 0.005, 0.005, 0.005, 0.005, 0.005]
EAN_RMSD_w = [0.18, 0.07, 0.05, 0.10, 0.05, 0.04, 0.03, 0.04, 0.03, 0.02, 0.01, 0.03, 0.02, 0.02, 0.01, 0.02, 0.06 ]
EAN_RMSD_s = [0.09, 0.01, 0.01, 0.02, 0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.005, 0.01, 0.005, 0.005, 0.005, 0.01, 0.02]

w = 1.0 # bar width
coast='open_sea'
nSub=len(OGS.P.basin_list)


for isub, sub in enumerate(OGS.P):
    print sub.name
    outfile=OUTFIG_DIR + "NRTvalidation_chlsup_" + sub.name + "_" + coast + ".png"
    #matplotlib.rc('xtick', labelsize=12)
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    # BIAS

    times, f0= F0.plotdata(F0.bias, sub.name,coast) ; ax1.bar(times,f0, width=w, bottom=0, align="center", color="grey",      label='1$^{st}$ day of forecast')
    times, f1= F1.plotdata(F1.bias, sub.name,coast);  ax1.bar(times,f1, width=w, bottom=0, align='center', color="skyblue",   label='2$^{nd}$ day of forecast')
    times, f2= F2.plotdata(F2.bias, sub.name,coast);  ax1.bar(times,f2, width=w, bottom=0, align='center', color="limegreen", label='3$^{rd}$ day of forecast')
    EAN_BIAS=np.zeros((len(times),nSub), dtype = np.float32)*np.nan
    for it, t in enumerate(times):
        if (( t.month <= 5 ) | ( t.month >= 10 )):
           EAN_BIAS[it,isub] = np.nan
        else:
           EAN_BIAS[it,isub] = EAN_BIAS_s[isub]
    ax1.plot(times,EAN_BIAS[:,isub],marker='_',color='red',linewidth=2.0, label="EAN")

    
    ax1.set_title(sub.extended_name + " (" + coast + ")",fontsize=14)
    ax1.set_ylabel('BIAS mg/m$^3$', fontsize=14)


    ax1.legend(bbox_to_anchor=(0.85,1), bbox_transform=pl.gcf().transFigure, fontsize = 11)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.grid(True)
    
    #RMS
    times, f0= F0.plotdata(F0.rmse, sub.name,coast) ; ax2.bar(times,f0, width=w, bottom=0, align="center", color="grey")
    times, f1= F1.plotdata(F1.rmse, sub.name,coast) ; ax2.bar(times,f1, width=w, bottom=0, align='center', color="skyblue")
    times, f2= F2.plotdata(F2.rmse, sub.name,coast) ; ax2.bar(times,f2, width=w, bottom=0, align='center', color="limegreen")

    EAN_RMSD=np.ndarray((len(times),nSub), dtype = np.float32)
    for it, t in enumerate(times):
        if (( t.month <= 5 ) | ( t.month >= 10 )):
           EAN_RMSD[it,isub] = EAN_RMSD_w[isub]
        else:
           EAN_RMSD[it,isub] = EAN_RMSD_s[isub]
    new_times, new_y = avoid_diagonals_in_jumps(times, EAN_RMSD[:,isub])
    ax2.plot(new_times, new_y,marker='_',color='red',linewidth=2.0)
    
    ax2.set_ylabel('RMS mg/m$^3$',fontsize=14)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.grid(True)


    # n. of valid points
    times, f0= F0.plotdata(F0.number, sub.name,coast) ; ax3.bar(times,f0, width=w, bottom=0, align="center", color="grey")
    times, f1= F1.plotdata(F1.number, sub.name,coast) ; ax3.bar(times,f1, width=w, bottom=0, align='center', color="skyblue")
    times, f2= F2.plotdata(F2.number, sub.name,coast) ; ax3.bar(times,f2, width=w, bottom=0, align='center', color="limegreen")
    

    ax3.set_ylabel('# of points',fontsize=14)
    ax3.xaxis.set_major_locator(mdates.MonthLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%b-%y"))

    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.savefig(outfile)
    pl.close(fig)

