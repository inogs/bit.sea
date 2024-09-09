import numpy as np
from bitsea.Sat import SatManager as Sat
import matplotlib
#matplotlib.use('Qt5Agg')
matplotlib.use('Agg')
import pylab as pl


Q10=np.load("/g100_scratch/userexternal/gbolzon0/V10C/bit.sea/Sat/KD490/KD490/2019/QI_threshold_1.0/rejected_counters.npy")
Q15=np.load("/g100_scratch/userexternal/gbolzon0/V10C/bit.sea/Sat/KD490/KD490/2019/QI_threshold_1.5/rejected_counters.npy")
Q20=np.load("/g100_scratch/userexternal/gbolzon0/V10C/bit.sea/Sat/KD490/KD490/2019/QI_threshold_2.0/rejected_counters.npy")
Q25=np.load("/g100_scratch/userexternal/gbolzon0/V10C/bit.sea/Sat/KD490/KD490/2019/QI_threshold_2.5/rejected_counters.npy")
Q30=np.load("/g100_scratch/userexternal/gbolzon0/V10C/bit.sea/Sat/KD490/KD490/2019/QI_threshold_3.0/rejected_counters.npy")


# nPoints_tot = Q10['tot']
# pl.close('all')
# fig,ax = pl.subplots()
#
#
# ax.plot(100*Q10['rejected']/nPoints_tot, 'r', label='1.0')
# ax.plot(100*Q15['rejected']/nPoints_tot, 'g', label='1.5')
# ax.plot(100*Q20['rejected']/nPoints_tot, 'b', label='2.0')
# ax.plot(100*Q25['rejected']/nPoints_tot, 'k', label='2.5')
# ax.plot(100*Q30['rejected']/nPoints_tot, 'm', label='3.0')
#
# ax.legend()
# ax.grid(True)
# #fig.show()






season = "winter"
filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/SAT/KD490/DT/DAILY/ORIG/20190301_cmems_obs-oc_med_bgc-transp_my_l3-multi-1km_P1D.nc"
#season = "summer"
#filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/SAT/KD490/DT/DAILY/ORIG/20190702_cmems_obs-oc_med_bgc-transp_my_l3-multi-1km_P1D.nc"


QI = Sat.readfromfile(filename,'QI_KD490')
VALUES_ORIG = Sat.readfromfile(filename,'KD490')


class sub_basin():
    def __init__(self,name,i_min,i_max, j_min, j_max, max_value_win, max_value_sum):
        self.name = name
        self.i_imin = i_min
        self.i_imax = i_max
        self.j_min = j_min
        self.j_max = j_max
        self.max_value_win = max_value_win
        self.max_value_sum = max_value_sum

S1 = sub_basin('ion2',1635,2167,1,666,0.06, 0.04)
S2 = sub_basin('nwm', 390, 1167,987,1480, 0.15, 0.06) #[-1,9, 40, 45 ]
S3 = sub_basin('lev', 2490, 3112, 198, 523, 0.06, 0.03 ) # [26,34,32, 35.3]



kd_min=0.021


for s in [ S1, S2, S3 ]:
    if season == "winter": max_value = s.max_value_win
    if season == "summer": max_value = s.max_value_sum
    
    bins=np.arange(0.016,max_value,0.001)
    VALUES = VALUES_ORIG.copy()
    LocalValues = VALUES[s.j_min:s.j_max, s.i_imin:s.i_imax]
    LocalQI =         QI[s.j_min:s.j_max, s.i_imin:s.i_imax]
    
    kd = LocalValues[LocalValues>0]
    Hist, bin_edges=np.histogram(kd,bins)
    
    #pl.close('all')
    #fig,ax=pl.subplots()
    #ax.bar(bin_edges[:-1],Hist, width=0.0005,align='edge',color='b')
    #fig.show()
    
    
    #low_values = (VALUES<=kd_min) & (VALUES>0)
    #VALUES[low_values] = kd_min
    THRESHOLD = 2.0
    bad = np.abs(QI) > THRESHOLD # 2.0
    #bad[low_values] = False
    
    VALUES[bad] = Sat.fillValue
    kd = LocalValues[LocalValues>0]
    Hist, bin_edges=np.histogram(kd,bins)
    
    
    
    
    print(kd.mean())
    fig,ax=pl.subplots()
    fig.set_size_inches(10,5)
    ax.bar(bin_edges[:-1],Hist, width=0.0005,align='edge', color='r', label="QI 2.0")
    ax.grid(True)
    
    
    THRESHOLD = 3.0
    VALUES = VALUES_ORIG.copy()
    LocalValues = VALUES[s.j_min:s.j_max, s.i_imin:s.i_imax]
    bad = np.abs(QI) > THRESHOLD # 2.0
    
    VALUES[bad] = Sat.fillValue
    kd = LocalValues[LocalValues>0]
    Hist, bin_edges=np.histogram(kd,bins)
    
    ax.bar(bin_edges[:-1],Hist, width=0.0003,align='edge',color='b', label="QI 3.0") 
    ax.set_title(s.name)
    ax.legend()
    ax.set_xlabel('Kd value')
    ax.set_ylabel('# accepted points')
    outfile = "kd_cmp_%s.%s.png" %(s.name, season)
    fig.savefig(outfile)


