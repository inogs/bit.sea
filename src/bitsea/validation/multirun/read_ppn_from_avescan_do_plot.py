import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates Figures IV.5 and others
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =False,
                                default = "none",
                                help = 'Input scanned ppn NetCDF files'
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = 'Fig',
                                help = 'Output images directory')


    return parser.parse_args()


args = argument()

import os
import scipy.io.netcdf as NC
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as pl
import math, operator
import glob
from dateutil.relativedelta import relativedelta
import datetime
from matplotlib.dates import date2num, num2date, YearLocator, MonthLocator, DateFormatter


# ATTENZIONE VA LANCIATO IL . ./config_MASKFILE.sh
#maskfile    = os.getenv("MASKFILE");
#ncIN=NC.netcdf_file(maskfile,"r")
#nav_lev = ncIN.variables['nav_lev'].data.copy()
#ncIN.close()

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')



DATADIR =args.inputdir
FIGDIR = args.outdir
os.system('mkdir -p ' + FIGDIR)

pl.close('all')
var='ppn'
SUBlist_subset=['alb','sww','swe','nwm','tyr','adn','ads','ion','lev','med']
subreg_mod=SUBlist_subset
#TEXTylabel="mgC/m2/day"
TEXTylabel="gC/m2/yr"

# define dates
dd=[]
#datelist=glob.glob(DATADIR + "ave." + in_yr + "*-12:00:00.col_integrals.nc")
datelist=glob.glob(DATADIR + "ave.*.col_integrals.nc")
datelist.sort()
for date in datelist:
    data=os.path.basename(date)
    dd.append(data[4:12])


nTimes=len(datelist)
TIMELIST = np.zeros(nTimes,np.float32)

iTime=0
DATE=[]
for data1 in dd:
    yr=data1[0:4]
    mn=data1[4:6]
    dy=data1[6:8]
    newdata=datetime.date(int(yr),int(mn),int(dy))
    datac=str(newdata)[2:4]+str(newdata)[5:7]+str(newdata)[8:10]
    TIMELIST[iTime] = newdata.toordinal() - datetime.datetime(1970,1,1).toordinal()
    DATE.append(data1)

Tlist=[]
for date in DATE:
  T=datetime.datetime.strptime(date,'%Y%m%d')
  Tlist.append(T)

# nc =  0 coast, 1 open_sea, 2 everywhere 
nc=2
# nd = 0 shallow (0-200), 1 deep (200-bottom)
nd=0

cf=0.365 #conversion factor for ppn: mgC m-2 d-1 = 365/1000 gC m-2 y-1

ppnM=np.zeros((len(subreg_mod),len(datelist)),np.float32)
ppnS=np.zeros((len(subreg_mod),len(datelist)),np.float32)


# READ data from REANALYSIS:
for imon,monfile in enumerate(datelist):
   P=NC.netcdf_file(monfile,"r")
   for ns, sub in enumerate(subreg_mod):
      ppnM[ns,imon]=P.variables['ppn'].data[ns,nc,nd,0] # leggo mean of ppn [mg/m2/d]
      ppnS[ns,imon]=P.variables['ppn'].data[ns,nc,nd,1] # leggo std

ppnMA=np.zeros((len(subreg_mod),len(datelist)/12),np.float32)
ppnSA=np.zeros((len(subreg_mod),len(datelist)/12),np.float32)
ppnMM=np.zeros((len(subreg_mod),12),np.float32)
ppnSM=np.zeros((len(subreg_mod),12),np.float32)


for ns, sub in  enumerate(subreg_mod):
# media annaule
  for ia in range (len(datelist)/12):
     ppnMA[ns,ia]=np.mean(ppnM[ns,ia*12:ia*12+12])*cf # annaul mean in gC/m2/y
     ppnSA[ns,ia]=np.mean(ppnS[ns,ia*12:ia*12+12])*cf
# climatologia mensile
  for im in range(0,12):
     ppnMM[ns,im]=np.mean(ppnM[ns,im::12]) # prendo uno ogni 12  MEDIE MENSILI in mg/m2/d
     ppnSM[ns,im]=np.mean(ppnS[ns,im::12])

# PLOT DEI ANNI IN gC/2/y
# build palette color of number(subreg_mod) colours
color=iter(cm.rainbow(np.linspace(0,1,len(subreg_mod)-2))) # escludo ADN e MED
Ianni=1999+np.arange(len(datelist)/12)
fig1=pl.figure(num=None,figsize=(8,6),dpi=100,facecolor='w',edgecolor='k')
ax1=fig1.add_subplot(111)
for ns,sub in enumerate(subreg_mod[:-1]):
  if ns <> 5 :# escludo ADN e MED
    print ns, sub
    c=next(color)
    pl.plot(Ianni,ppnMA[ns,:],'-',c=c,label=sub)
ax1.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
ax1.grid(color='k',linestyle='--')    
ax1.ticklabel_format(fontsize=12)
ax1.set_ylabel('[PPN gC/m2/y]',fontsize=12)
pl.rc('xtick', labelsize=12)
pl.rc('ytick', labelsize=12)
leg = pl.gca().get_legend()
ltext  = leg.get_texts()
pl.setp(ltext,fontsize=12)
pl.show(block=False)
nomefig= FIGDIR + "Plot_0_200vert_integ_ppn_offshore_annual_timeserie.png"
fig1.savefig(nomefig, dpi=fig1.dpi)

# PLOT DEI MESI in mg/m2/d
x_ticks= []
x_labels=[]
color=iter(cm.rainbow(np.linspace(0,1,len(subreg_mod)-2)))
Imesi=['J','F','M','A','M','J','J','A','S','O','N','D']
x=np.arange(12)
for im,lab in enumerate(Imesi):
 x_labels.append(lab)
 x_ticks.append(im)

fig1=pl.figure(num=None,figsize=(8,6),dpi=100,facecolor='w',edgecolor='k') # escludo ADN e MED
ax1=fig1.add_subplot(111)
for ns,sub in enumerate(subreg_mod[:-1]):
   if ns <> 5:  # escludo ADN e MED
    print ns, sub
    c=next(color)
    pl.plot(x,ppnMM[ns,:],'-',c=c,label=sub)
ax1.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
ax1.grid(color='k',linestyle='--')
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_labels)
ax1.ticklabel_format(fontsize=12)
ax1.set_ylabel('[PPN mgC/m2/d]',fontsize=12)
pl.rc('xtick', labelsize=10)
pl.rc('ytick', labelsize=10)
leg = pl.gca().get_legend()
ltext  = leg.get_texts()
pl.setp(ltext,fontsize=12)
pl.show(block=False)

nomefig= FIGDIR + "Plot_0_200vert_integ_ppn_offshore_monthly_climatol.png"
fig1.savefig(nomefig, dpi=fig1.dpi)



# plot





