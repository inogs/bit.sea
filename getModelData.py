import Float_Manager
import numpy as np
import postproc
import postproc.genUserDateList as DL
from postproc.Timelist import *
from postproc.requestors import Daily_req

def dump_punti_for_aveScan(Floatlist,filename):
    LINES=[]
    LINES.append('NOME    Longitutine N    Latitudine E\n')
    for f in Floatlist:
        line = "%s %g %g \n" %(f.wmo,f.lat, f.lon)
        LINES.append(line)
    
    F = open(filename, "w")
    F.writelines(LINES)
    F.close()

INPUTDIR="/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/FAKE/AVE/"
MODEL_PROFILES_DIR="/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/FAKE/PROFILES/"
PUNTI_DIR="/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/FAKE/PUNTI/"
TMPSDIR = "/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/FAKE/"
# ModTimes = DL.getTimeList('20150601-12:00:00', '20150620-12:00:00',"days = 1")
# for t in ModTimes[:0]:
#     filename=INPUTDIR + t.strftime("ave.%Y%m%d-%H:%M:%S.nc")
#     F=open(filename,'w')
#     F.close()

TL = TimeList('20150601-12:00:00', '20150620-12:00:00', INPUTDIR,"ave*.nc",'postproc/IOnames.xml')

var = 'NITRATE'
TI = Float_Manager.Time_Interval('20150520','20150830','%Y%m%d')
R = Float_Manager.Region(-6,36,30,46)
FLOAT_LIST=Float_Manager.FloatSelector(None, TI, R)

Coupled_List=[]
for ir, req in enumerate(TL.getOwnList()):
    LIST_of_REQ=[]
    for f in FLOAT_LIST:
        if (f.time >= req.starttime) & (f.time <= req.endtime) :
            LIST_of_REQ.append(f)
    if (len(LIST_of_REQ) >0 ): Coupled_List.append((TL.Timelist[ir],LIST_of_REQ))
    


JOB_LINES=[]
for it, t in enumerate(Coupled_List):
    Model_time        = t[0]
    INTERESTED_FLOATS = t[1]
    
    outpuntifile= PUNTI_DIR + "punti_" + Model_time.strftime("%Y%m%d") + ".dat" #punti_20150416.dat
    dump_punti_for_aveScan(INTERESTED_FLOATS, outpuntifile)
    line = 'python aveScan.py '   + \
        ' -i '  + INPUTDIR  +  \
        ' -t '  +  TMPSDIR  + \
        ' -o '  + MODEL_PROFILES_DIR  + \
        ' -d VarDescriptorB.xml ' + \
        ' -p ' + outpuntifile + '\n'
    JOB_LINES.append(line)
    F=file('jobProfiler.pbs','w')
    F.writelines(JOB_LINES)
    F.close()

# Provare lo stesso per le fisiche
