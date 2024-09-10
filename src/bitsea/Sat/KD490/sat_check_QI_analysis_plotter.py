import pickle
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl
from bitsea.postproc import masks
from datetime import datetime
from bitsea.layer_integral import coastline

class container():
    def __init__(self, I,J,values, QI_old, QI_new):
        self.I = I
        self.J = J
        self.values = values
        self.QI_old = QI_old
        self.QI_new = QI_new




pl.close('all')
INPUTDIR="CHL/2019/QI_threshold_2"

if True:
    A=np.load(INPUTDIR + '/rejected_counters.npy')
    fig,ax=pl.subplots()

    ax.plot(A['only_clim_reject'],'r',label='rejected CLIM')
    ax.plot(A['only_qi_reject'],'b',label='rejected QI')
    ax.legend()
    fig.show()


    fig,ax=pl.subplots()
    nPoints_tot = A['tot']
    ax.plot(100*A['only_clim_reject']/nPoints_tot,'r',  label='% rejected only CLIM')
    ax.plot(100*A['only_qi_reject']/nPoints_tot,  'b',  label='% rejected only QI')
    ax.plot(100*A['both_reject']/nPoints_tot,  'g',     label='% common rejection')
    ax.legend()
    fig.show()

    ratio = A['only_qi_reject'].sum()/A['only_clim_reject'].sum()
    print ("QI_rejected/CLIM_rejected=",ratio)




coastlinelon, coastlineclat = coastline.get()
maskSat = getattr(masks,'SAT1km_mesh')




fid = open(INPUTDIR + '/rejected_only_old.pkl','rb')
DAILY_REJECT_ONLY_OLD=pickle.load(fid)
fid.close()

fid = open(INPUTDIR + '/rejected_only_new.pkl','rb')
DAILY_REJECT_ONLY_NEW = pickle.load(fid)
fid.close()



datestart    = '20190701'
d            = datetime.strptime(datestart,"%Y%m%d")
iFrame_start = int(d.strftime('%j')) -1
iFrame_end   = iFrame_start + 1

I_clim = np.array([],int)
J_clim = np.array([],int)
I_QI   = np.array([],int)
J_QI   = np.array([],int)


VALUES = np.array([],np.float32)
for C in DAILY_REJECT_ONLY_OLD[iFrame_start: iFrame_end]:
    I_clim = np.concatenate((I_clim,C.I))
    J_clim = np.concatenate((J_clim,C.J))
    VALUES = np.concatenate((VALUES, C.values))

for C in DAILY_REJECT_ONLY_NEW[iFrame_start: iFrame_end]:
    I_QI = np.concatenate((I_QI,C.I))
    J_QI = np.concatenate((J_QI,C.J))



fig,ax = pl.subplots()
ax.plot(coastlinelon, coastlineclat,'k')
ax.plot(maskSat.lon[I_clim],maskSat.lat[J_clim],'r.', markersize=4, label='rejected only CLIM')
ax.plot(maskSat.lon[I_QI],maskSat.lat[J_QI],'b.',markersize=4, label='rejected only QI')
ax.legend()

fig.show()



hist, bin_edges=np.histogram(VALUES,bins=np.arange(0,1,0.02))
for i in hist[:10]: print ("%f %% " %(100*i/len(VALUES)))

#C=DAILY_REJECT_ONLY_NEW[iFrame_start]
#fig,ax=pl.subplots()
#ax.plot(C.QI_old,'.')
#fig.show()

# search reason of CLIM rejection in  DAILY_REF_MEAN and DAILY_REF_STD
# ii=(VALUES > 0.34) & (VALUES <=0.36)
# Itop=I_clim[ii]
# Jtop=J_clim[ii]
# nTop=4593
#
# MEAN_CLIM = np.zeros((nTop),np.float32)
# MEAN_STD  = np.zeros((nTop),np.float32)
# for k in range(nTop):
#     j = Jtop[k]
#     i = Itop[k]
#     MEAN_CLIM[k] = DAILY_REF_MEAN[j,i]
#     MEAN_STD[k] =  DAILY_REF_STD[j,i]
#
# fig,ax=pl.subplots()
# ax.plot(MEAN_CLIM,'.')
# fig.show()
# fig,ax=pl.subplots()
# ax.plot(MEAN_STD,'.')
# fig.show()
# -----------------------------------------------


