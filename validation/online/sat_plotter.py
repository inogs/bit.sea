import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates time series png files for mdeeaf web site.
    Each file has 3 subplots: bias, rms and number of points.
    Only the daily matchup sat/model is taken in account
    The time window is 18 months.
    See https://matplotlib.org/3.1.0/gallery/color/named_colors.html for colors
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--date','-d',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'input directory with all the Validation_*nc, operationally inpdir/SAT_VALIDATION/')

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
import matplotlib.pyplot as pl
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
INPUT_DIR        = addsep(args.inputdir)
OUTFIG_DIR       = addsep(args.outdir)



prefix="Validation_f1_20190723_on_daily_Sat."
F0 = timelistcontainer(TI,INPUT_DIR, 'Validation_f1_*on_daily_Sat*', prefix=prefix)
F1 = timelistcontainer(TI,INPUT_DIR, 'Validation_f2_*on_daily_Sat*', prefix=prefix)
F2 = timelistcontainer(TI,INPUT_DIR, 'Validation_f3_*on_daily_Sat*', prefix=prefix)



#alb swm1 swm2 nwm tyr1 tyr2 adr1 adr2 aeg ion1 ion2 ion3 lev1 lev2 lev3 lev4 med

EAN_BIAS_s = [0.008, 0.005, 0.005, 0.005, 0.005, 0.005, -0.01, -0.007, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]

EAN_RMSD_w = [0.10, 0.05, 0.05, 0.05, 0.03, 0.03, 0.04, 0.04, 0.03, 0.02, 0.02, 0.03, 0.02, 0.02, 0.01, 0.01, 0.04]
EAN_RMSD_s = [0.05, 0.01, 0.006, 0.008, 0.007, 0.005, 0.01, 0.01, 0.008, 0.005, 0.005, 0.007, 0.005, 0.005, 0.005, 0.008, 0.01]


w = 0.9 # bar width
coast='open_sea'
nSub=len(OGS.P.basin_list)


for isub, sub in enumerate(OGS.P):
    if sub.name == 'atl': continue
    print sub.name
    outfile=OUTFIG_DIR + "NRTvalidation_chlsup_" + sub.name + "_" + coast + ".png"
    #matplotlib.rc('xtick', labelsize=12)
    fig, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, figsize=(15,10)) 
    # BIAS

    times, f0= F0.plotdata(F0.bias, sub.name,coast) ; ax1.plot(times,f0,'-k' , label='1-day lead time fc') # label='1$^{st}$ day of forecast')
    times1=times
    times, f1= F1.plotdata(F1.bias, sub.name,coast);  ax1.plot(times,f1,'--b', label='2-day lead time fc') # label='2$^{nd}$ day of forecast')  
    times, f2= F2.plotdata(F2.bias, sub.name,coast);  ax1.plot(times,f2,':g' , label='3-day lead time fc') # label='3$^{rd}$ day of forecast')
    times2=times
    times1.extend(times2[5:])

    EAN_BIAS=np.zeros((len(times1),nSub), dtype = np.float32)*np.nan
    for it, t in enumerate(times1):
        if (( t.month <= 5 ) | ( t.month >= 10 )):
           EAN_BIAS[it,isub] = np.nan
        else:
           EAN_BIAS[it,isub] = EAN_BIAS_s[isub]
    ax1.plot(times1,EAN_BIAS[:,isub],marker='_',color='red',linewidth=2.0, label="EAN")

    
    ax1.set_title(sub.extended_name + " (" + coast + ")",fontsize=14)
    ax1.set_ylabel('BIAS mg/m$^3$', fontsize=14)


    ax1.legend(bbox_to_anchor=(1.005,.95), loc="lower right" , fontsize = 11)#bbox_transform=pl.gcf().transFigure,
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.grid(True)
    
    #RMS
    F0 = timelistcontainer(TI,INPUT_DIR, 'Validation_f1_*on_daily_Sat*', prefix=prefix)
    times, f0= F0.plotdata(F0.rmse, sub.name,coast) ; 
    ax2.plot(times,f0,'-k')
    times, f1= F1.plotdata(F1.rmse, sub.name,coast) ; ax2.plot(times,f1,'--b') 
    times, f2= F2.plotdata(F2.rmse, sub.name,coast) ; ax2.plot(times,f2,':g')

    EAN_RMSD=np.ndarray((len(times1),nSub), dtype = np.float32)
    for it, t in enumerate(times1):
        if (( t.month <= 5 ) | ( t.month >= 10 )):
           EAN_RMSD[it,isub] = EAN_RMSD_w[isub]
        else:
           EAN_RMSD[it,isub] = EAN_RMSD_s[isub]
    new_times, new_y = avoid_diagonals_in_jumps(times1, EAN_RMSD[:,isub])
    ax2.plot(new_times, new_y,marker='_',color='red',linewidth=2.0)
    
    ax2.set_ylabel('RMS mg/m$^3$',fontsize=14)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.grid(True)


    # n. of valid points
    times, f0= F0.plotdata(F0.number, sub.name,coast) ; ax3.bar(times,f0, width=w, bottom=0, align="center", color="grey")
    

    ax3.set_ylabel('# of points',fontsize=14)
    ax3.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=0, interval=2))
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%d-%b-%y"))

    fig.subplots_adjust(hspace=0.1)
    pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.savefig(outfile)
    pl.close(fig)

