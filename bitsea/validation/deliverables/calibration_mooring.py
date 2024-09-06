from instruments import mooring
from instruments.var_conversions import MOORINGVARS
from commons.time_interval import TimeInterval
import numpy as np
from profiler import *
from basins.region import Rectangle
from commons.layer import Layer
import matplotlib.pyplot as pl

import os
import scipy.io.netcdf as NC
M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)

maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

modelvarname = 'O2o' 
region = Rectangle(-6,36,30,46)

rE1M3A = Rectangle(25,26,35,36)
rPylos = Rectangle(21,22,36,37)
rCor   = Rectangle(9,10,43,44)
rNADR  = Rectangle(12,13,44,45)


L20 = Layer(20,22)
L50 = Layer(45,55)
L75 = Layer(70,80)
L100 =Layer(95,105)
MOORING_LIST =['E1M3A','Pylos','COR','NADR']
REGION__LIST =[rE1M3A, rPylos, rCor, rNADR]

MOORING_LIST =['Pylos']
REGION__LIST =[ rPylos ]

#MOORING_LIST =['E1M3A']
#REGION__LIST =[ rE1M3A ]

MOORING_LIST =['COR']
REGION__LIST =[ rCor ]

MOORING_LIST =['NADR']
REGION__LIST =[rNADR]


STRINGLIST=['20m','50m','75m','100m']


var = MOORINGVARS[modelvarname]
ProfileList=mooring.MooringSelector(var, T_INT, region)
# nP = len(ProfileList)
# Lon=np.zeros((nP,),np.float32)
# Lat=np.zeros((nP,),np.float32)
# for ip,p in enumerate(ProfileList):
#     Lon[ip] = p.lon
#     Lat[ip] = p.lat

#
for imooring, mooringstring in enumerate(MOORING_LIST):
    print mooringstring
#    ProfileList=mooring.MooringSelector('CPHL', T_INT, REGION__LIST[imooring])
    ProfileList=mooring.MooringSelector(var, T_INT, REGION__LIST[imooring])
    m = M.getMatchups(ProfileList, nav_lev, modelvarname, read_adjusted=True)
#    fig, ax = m.densityplot(30)
#    ax.set_ylabel('CHL Ref')
#    ax.set_xlabel('CHL Model')
#    filename = mooringstring + "_density.png"
#    fig.savefig(filename)
#    pl.clf()

#    import sys
#    sys.exit() 
    fig1=pl.figure(num=None,facecolor='w',edgecolor='k')
    for ilayer, layer in enumerate([L20, L50, L75, L100]): 
        mlayer = m.subset(layer)
        ax1=fig1.add_subplot(2,2,ilayer+1)
        pl.plot(mlayer.Time,mlayer.Ref,'.b')
        pl.plot(mlayer.Time,mlayer.Model,'.r')
#        pl.ylim(0,0.4)
        pl.title(layer.string())
#        ax1.set_ylabel('CHL [mg/m^3]')
        ax1.set_ylabel('OXY [mmol/m3]')
        ax1.set_xticklabels("")

#        ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,4,7,10],bymonthday=10))
#        ax1.xaxis.set_major_locator(mdates.MonthLocator()
#        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b%y"))
    pl.show(block=False) 
    filename = "plot_OXY_" + mooringstring + ".png"
    fig1.savefig(filename)
#    pl.clf()


