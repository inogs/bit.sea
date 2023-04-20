from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from postproc import masks
import numpy as np
import os

from Sat import SatManager as Sat

ORIGDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C//SAT/KD490/DT/DAILY/ORIG/" 
CLIM_FILE = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/KD490_Climatology_1km.nc"
maskSat = getattr(masks,'SAT1km_mesh')

class container():
    def __init__(self, I,J,values, QI_old, QI_new):
        self.I = I
        self.J = J
        self.values = values
        self.QI_old = QI_old
        self.QI_new = QI_new


reset = False

Timestart="20200101"
Time__end="20210101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')

nFrames = TL_orig.nTimes
REJECT=np.zeros((nFrames,),dtype=[('tot',int),('both_reject',int),('only_clim_reject',int), ('only_qi_reject',int)])
DAILY_REJECT_ONLY_OLD = []
DAILY_REJECT_ONLY_NEW = []

MEAN,STD = Sat.readClimatology(CLIM_FILE)


for iTime, filename in enumerate(TL_orig.filelist):
    print(iTime)

    julian = int( TL_orig.Timelist[iTime].strftime("%j") )
    if julian==366: julian=365 
    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    

    VALUES_IN = Sat.readfromfile(filename,'KD490')
    QI        = Sat.readfromfile(filename,'QI_KD490')
    cloudsLandTIME = VALUES_IN         == Sat.fillValue
    cloudlandsCLIM = DAILY_REF_MEAN == Sat.fillValue
    
    CHL_OUT = VALUES_IN.copy()
    CHL_OUT[cloudsLandTIME] = Sat.fillValue
    CHL_OUT[cloudlandsCLIM] = Sat.fillValue
    counter_refNAN = (~cloudsLandTIME & cloudlandsCLIM).sum(axis=None)
    
    
    outOfRange = np.abs(VALUES_IN - DAILY_REF_MEAN) > DAILY_REF_STD *2.0
    outOfRange[cloudsLandTIME | cloudlandsCLIM ] = False
    Elim = outOfRange | (~cloudsLandTIME & cloudlandsCLIM)
    
    
    
    CHL_OUT[outOfRange] = Sat.fillValue 


    qi_cloud = VALUES_IN == Sat.fillValue
    ############################################################
    outOfRange_qi = np.abs(QI) > 2
    ############################################################
    outOfRange_qi[qi_cloud] = False
    

    
    

    out_both     =  Elim &  outOfRange_qi
    out_only_old =  Elim & ~outOfRange_qi
    out_only_new = ~Elim &  outOfRange_qi

    counter_qi_out = outOfRange_qi.sum(axis=None)
    
    values_only_old = VALUES_IN[out_only_old]
    values_only_new = VALUES_IN[out_only_new]

    REJECT[iTime]['tot']= (~cloudsLandTIME).sum()
    REJECT[iTime]['both_reject'] = out_both.sum()
    REJECT[iTime]['only_clim_reject']   = out_only_old.sum()
    REJECT[iTime]['only_qi_reject']     = out_only_new.sum()
    
    J,I = np.nonzero(out_only_old)
    QI_old = np.zeros_like(I)
    daily_only_old_object = container(I,J,VALUES_IN[out_only_old], QI_old, QI[out_only_old])
    DAILY_REJECT_ONLY_OLD.append(daily_only_old_object)
    
    J,I = np.nonzero(out_only_new)
    QI_old = (VALUES_IN[out_only_new]-DAILY_REF_MEAN[out_only_new])/DAILY_REF_STD[out_only_new]
    daily_only_new_object  = container(I,J,VALUES_IN[out_only_new],QI_old,QI[out_only_new])
    DAILY_REJECT_ONLY_NEW.append(daily_only_new_object)
    
    #import sys
    #sys.exit()
    
import pickle
fid = open('rejected_only_old.pkl','wb')
pickle.dump(DAILY_REJECT_ONLY_OLD, fid)
fid.close()

fid = open('rejected_only_new.pkl','wb')
pickle.dump(DAILY_REJECT_ONLY_NEW, fid)
fid.close()

np.save('rejected_counters',REJECT)

