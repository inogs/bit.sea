import numpy as np
import argparse
import matplotlib as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as pl
from basins import V2 as OGS
from commons import timerequestors
from timeseries.plot import read_pickle_file
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import writetable
from datetime import datetime


def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--indir', '-i',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Input dir'''
                            )
    parser.add_argument(   '--reandir', '-r',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Reanalysis'''
                            )
    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )




    return parser.parse_args()


args = argument()

#REANDIR = open(args.reandir)
REANDIR  = args.reandir
INPUTDIR = args.indir
OUTDIR   = args.outdir

#TI=TimeInterval("19990101","20180101",'%Y%m%d')
#TL=TimeList.fromfilenames(TI, REANDIR, "ave*nc")
#ML=TL.getMonthlist()
#MClim=timerequestors.Clim_month

istart=14 #index from JAN 2017

VARLIST=[ 'pCO2','CO2airflux']
UNITS  =[ '$ [ \mu  atm] $ ' , '$ [ mmol \, m^{-2} \, day^{-1} ] $']
StatDescr   = ["Mean", "Std", "min", "p05", "p25", "p50", "p75", "p95", "max"] # 0=MEAN, 2=MIN, 8=MAX
COASTNESS_LIST=['coast', 'open_sea','everywhere'] #1='open_sea'
nSUB = len(OGS.P.basin_list)
#CLIM_MONTHLY = np.zeros((nSUB, 12),np.float32)*np.nan
REAN_MEAN = np.zeros((nSUB, 12),np.float32)*np.nan
REAN_MIN  = np.zeros((nSUB, 12),np.float32)*np.nan
REAN_MAX  = np.zeros((nSUB, 12),np.float32)*np.nan
MOD_MEAN_month = np.zeros((nSUB, 12),np.float32)*np.nan

#MOD_MEAN = np.zeros((nSUB, TL.nTimes),np.float32)*np.nan

########## MONTHS DEFINITION ##########
Imesi=['J','F','M','A','M','J','J','A','S','O','N','D']
x=np.arange(12)
x_ticks= []
x_labels=[]
for im,lab in enumerate(Imesi):
 x_labels.append(lab)
 x_ticks.append(im)

x_months=[datetime(2017,m+1,15) for m in range(12)]
########################################

for ivar, var in enumerate(VARLIST):
    print var
    REAN_MEAN = np.zeros((nSUB, 12),np.float32)*np.nan
    REAN_MIN  = np.zeros((nSUB, 12),np.float32)*np.nan
    REAN_MAX  = np.zeros((nSUB, 12),np.float32)*np.nan

    filename_REAN=REANDIR + var + ".pkl"
    REAN , TL_rean = read_pickle_file(filename_REAN)
    
    filename=INPUTDIR + var + ".pkl"
    data, TL = read_pickle_file(filename)
# from 2017 --> istart=14
    TL.Timelist=TL.Timelist[istart:]
    TL.filelist=TL.filelist[istart:]
    TL.nTimes=len(TL.filelist)
    data=data[istart:,:,:,:,:]

    MOD_MEAN = np.zeros((nSUB, TL.nTimes),np.float32)*np.nan

#    TL_rean.getMonthlist() #list of monthly req

    for isub, sub in enumerate(OGS.P):
        print sub.name
        MOD_MEAN[isub,:]=data[:,isub,1,0,0].mean()
        for imonth in range(12):
            MClim_req=timerequestors.Clim_month(imonth+1)  # for REAN
            ii,w=TL_rean.select(MClim_req)
        
            req=timerequestors.Monthly_req(2017,imonth+1)  # for MODEL
            im,wm = TL.select(req)
        
            REAN_MEAN[isub,imonth]=REAN[ii,isub,1,0,0].mean()
            REAN_MAX[isub,imonth]=REAN[ii,isub,1,0,8].max()
            REAN_MIN[isub,imonth]=REAN[ii,isub,1,0,2].min()
            MOD_MEAN_month[isub,imonth]=data[im,isub,1,0,0].mean() 

#    rows_names_list=[sub.name for sub in OGS.P]
#    column_names_list=[str(i) for i in range(1,13)]
#    outfile_mean="CLIM_REAN_MEAN_%s.txt" %(var) 
#    outfile_min="CLIM_REAN_MIN_%s.txt" %(var)
#    outfile_max="CLIM_REAN_MAX_%s.txt" %(var)
#    writetable(outfile_mean, REAN_MEAN, rows_names_list, column_names_list)
#    writetable(outfile_min, REAN_MAX, rows_names_list, column_names_list)
#    writetable(outfile_max, REAN_MIN, rows_names_list, column_names_list)
    
#    fig, axs = pl.subplots(2,2, facecolor='w', edgecolor='k')
#    for isub, sub in enumerate(OGS.P):
        fig, axs = pl.subplots()
        axs.plot(x_months,REAN_MEAN[isub,:],color='r',marker="D",label="mean (min and max) of REAN MED-BIO-006-008")
        axs.plot(x_months,REAN_MAX[isub,:],color='r',linestyle='--')
        axs.plot(x_months,REAN_MIN[isub,:],color='r',linestyle='--')

#        axs.plot(x,MOD_MEAN[isub,:],color='k',marker="o",label="mod year 2017")
#        axs.plot(TL.Timelist[istart:],data[istart:,isub,1,0,0],color='k',marker="o",label="mod year 2017")
        axs.plot(TL.Timelist,data[:,isub,1,0,0],color='k',marker=".",label="mod year 2017")
        axs.plot(x_months,MOD_MEAN_month[isub,:],color='k',linestyle='--')

        axs.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
#        axs.set_xticks(x_ticks)
#        axs.set_xticklabels(x_labels)
#        axs.ticklabel_format(fontsize=12)
        ylab = var + " " +  UNITS[ivar]
        axs.set_ylabel(ylab)
        axs.set_title(sub.name )
        axs.grid(color='k',linestyle='--')
        outfile = OUTDIR + var + "_" + sub.name + '_REAN_monthly_tseries_Fig4.20.png'
        fig.savefig(outfile, dpi=150)

import sys
sys.exit()

socat=np.loadtxt("monthly_clim_socat.txt",skiprows=1,usecols=range(1,13))
model=np.loadtxt("monthly_pCO2_TEST_04.txt",skiprows=1,usecols=range(1,13))

nSUB = len(OGS.P.basin_list)

x_ticks= []
x_labels=[]
color=iter(cm.rainbow(np.linspace(0,1,nSUB-1)))
color=iter(cm.rainbow(np.linspace(0,1,4)))
K=1.054 # Fattore correttivo di attualizzazione per la pCO2

Imesi=['J','F','M','A','M','J','J','A','S','O','N','D']
x=np.arange(12)
for im,lab in enumerate(Imesi):
 x_labels.append(lab)
 x_ticks.append(im)

fig, axs = pl.subplots(2,2, facecolor='w', edgecolor='k')
axs = axs.ravel()

for ns,sub in enumerate(OGS.P.basin_list[:4]):
            print ns, sub
            c=next(color)
            axs[0].plot(x,model[ns,:],'-',c=c,label=sub.name,linewidth=2.0)
            axs[0].plot(x,socat[ns,:]*K,'--',c=c,linewidth=2.0)
            axs[0].legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
            axs[0].grid(color='k',linestyle='--')
            axs[0].set_xticks(x_ticks)
            axs[0].set_xticklabels(x_labels)
            axs[0].ticklabel_format(fontsize=12)
            axs[0].set_ylabel("$\mu  atm$")

color=iter(cm.rainbow(np.linspace(0,1,4)))

for ns,sub in enumerate(OGS.P.basin_list[4:8]):
            js=ns+4
            print ns, sub
            c=next(color)
            axs[1].plot(x,model[js,:],'-',c=c,label=sub.name,linewidth=2.0)
            axs[1].plot(x,K*socat[js,:],'--',c=c,linewidth=2.0)
            axs[1].legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
            axs[1].grid(color='k',linestyle='--')
            axs[1].set_xticks(x_ticks)
            axs[1].set_xticklabels(x_labels)
            axs[1].ticklabel_format(fontsize=12)

