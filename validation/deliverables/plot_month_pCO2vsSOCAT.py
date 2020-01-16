import numpy as np
import matplotlib as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as pl
from basins import V2 as OGS

socat=np.loadtxt("monthly_clim_socat.txt",skiprows=1,usecols=range(1,13))
model=np.loadtxt("monthly_2017_surf/monthly_pCO2.txt",skiprows=1,usecols=range(1,13))

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
            axs[0].plot(x,socat[ns,:]*K,'--',c=c,linewidth=2.0) #,marker="o")
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
            axs[1].plot(x,K*socat[js,:],'--',c=c,linewidth=2.0) #,marker="o")
            axs[1].legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
            axs[1].grid(color='k',linestyle='--')
            axs[1].set_xticks(x_ticks)
            axs[1].set_xticklabels(x_labels)
            axs[1].ticklabel_format(fontsize=12)

color=iter(cm.rainbow(np.linspace(0,1,4)))

for ns,sub in enumerate(OGS.P.basin_list[8:12]):
            js=ns+8
            print ns, sub
            c=next(color)
            axs[2].plot(x,model[js,:],'-',c=c,label=sub.name,linewidth=2)
            axs[2].plot(x,K*socat[js,:],'--',c=c,linewidth=2) #,marker="o")
            axs[2].legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
            axs[2].grid(color='k',linestyle='--')
            axs[2].set_xticks(x_ticks)
            axs[2].set_xticklabels(x_labels)
            axs[2].ticklabel_format(fontsize=12) 
            axs[2].set_ylabel("$\mu  atm$")

color=iter(cm.rainbow(np.linspace(0,1,4)))
for ns,sub in enumerate(OGS.P.basin_list[12:16]):
            js=ns+12
            print ns, sub
            c=next(color)
            axs[3].plot(x,model[js,:],'-',c=c,label=sub.name,linewidth=2)
            axs[3].plot(x,K*socat[js,:],'--',c=c,linewidth=2) #,marker="o")
            axs[3].legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
            axs[3].grid(color='k',linestyle='--')
            axs[3].set_xticks(x_ticks)
            axs[3].set_xticklabels(x_labels)
            axs[3].ticklabel_format(fontsize=12)
#            axs[3].set_ylabel("$\mu  atm$")


fig.suptitle('BFMv5 pCO2 (solid) vs SOCAT pCO2 (dots)')
fig.savefig('pCO2_monthly_tseries_Fig4.20.png', dpi=150)

import sys
sys.exit()

fig1, ax1 = pl.subplots(2, 2)
fig1=pl.figure(num=None,figsize=(8,6),dpi=100,facecolor='w',edgecolor='k')
for i in [0,1]:
    for j in [0,1]:
        color=iter(cm.rainbow(np.linspace(0,1,nSUB-1)))
        for ns,sub in enumerate(OGS.P.basin_list[:-1]):
            print ns, sub
            c=next(color)
            ax1[i,j].plot(x,model[ns,:],'-',c=c,label=sub)
#ax1=fig1.add_subplot(111)
#for ns,sub in enumerate(OGS.P.basin_list[:-1]):
#    print ns, sub
#    c=next(color)
#    pl.plot(x,model[ns,:],'-',c=c,label=sub)
#    pl.plot(x,socat[ns,:],':',c=c)

ax1.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
ax1.grid(color='k',linestyle='--')
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_labels)
ax1.ticklabel_format(fontsize=12)
fig1.show()
