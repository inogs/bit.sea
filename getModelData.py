import Float_Manager
import numpy as np
def Time_Coupler(MeasTimes, ModelTimes):
    Merged_Timelist=[]
    Meas_Index = np.zeros(10,np.float32)
    Model_Index = np.zeros(10,np.float32)
    return Merged_Timelist,Meas_Index,Model_Index


def dump_punti_for_aveScan(Floatlist,filename):
    return 1


var = 'NITRATE'
TI = Float_Manager.Time_Interval('20150520','20150830','%Y%m%d')
R = Float_Manager.Region(-6,36,30,46)



FLOAT_LIST=Float_Manager.FloatSelector(var=None, TI, R)
nFloats = len(FLOAT_LIST)
Float_Times =[]
for f in FLOAT_LIST: Float_Times.append(f.time)

Model_Times = []
MODEL_PROFILES_DIR='/some/path'

Merged_List, Float_Index, Model_Index= Time_Coupler(Float_Times, Model_Times)

JOB_LINES=[]
for it, t in enumerate(Merged_List):
    INTERESTED_FLOATS=[]
    for i in np.arange(nFloats):
        if Float_Index[i] == it: INTERESTED_FLOATS.append(FLOAT_LIST[i])
    
    outpuntifile="punti_" + t.strftime("%Y%m%d") + ".dat" #punti_20150416.dat
    dump_punti_for_aveScan(INTERESTED_FLOATS, outpuntifile)
    line = 'python aveScan.py -o ' + MODEL_PROFILES_DIR  + ' -p ' + outpuntifile
    JOB_LINES.append(line)


Forcing_Times=[]    
Merged_List, Float_Index, Model_Index= Time_Coupler(Float_Times, Forcing_Times)