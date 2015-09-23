import Float_Manager
import numpy as np
import postproc
from postproc.Timelist import *
from postproc.requestors import Daily_req
from region import *
from matchup import *
import scipy.io.netcdf as NC
import os,sys
from shared_data import *


def dump_punti_for_aveScan(Floatlist,filename):
    LINES=[]
    LINES.append('NOME    Longitutine E    Latitudine N \n')
    for f in Floatlist:
        line = "%s %g %g \n" %(f.wmo,f.lon, f.lat)
        LINES.append(line)
    
    F = open(filename, "w")
    F.writelines(LINES)
    F.close()

def write_files_for_profiling(filename, Coupled_List,floatlist):
    JOB_LINES=[]
    
    JOB_LINES.append("export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc \n")
    JOB_LINES.append("export SUBMASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/POSTPROC/submask.nc \n")
    JOB_LINES.append("cd postproc \n")
    for t in Coupled_List:
        Model_time        = t[0]
        INTERESTED_FLOATS = [floatlist[k] for k in t[1]] #t[1]
        
        outpuntifile= PUNTI_DIR + "punti_" + Model_time.strftime("%Y%m%d") + ".dat" #punti_20150416.dat
        dump_punti_for_aveScan(INTERESTED_FLOATS, outpuntifile)
        line = 'python aveScan.py '   + \
            ' -l '  + Model_time.strftime("ave.%Y%m%d*")  + \
            ' -i '  + INPUTDIR  +  \
            ' -t '  + TMPSDIR  + \
            ' -o '  + BASEDIR  + \
            ' -d VarDescriptorB.xml ' + \
            ' -p ' + outpuntifile + '\n'
        JOB_LINES.append(line)
    
    F=file(filename,'w')
    F.writelines(JOB_LINES)
    F.close()



for DIR in [MODEL_PROFILES_DIR,TMPSDIR, PUNTI_DIR]: 
    os.system("mkdir -p " + DIR)

#  -----------------  init ----------------------
TL = TimeList(DATESTART, DATE__END, INPUTDIR,"ave*.nc",'postproc/IOnames.xml')
TI = Float_Manager.Time_Interval(DATESTART,DATE__END,'%Y%m%d-%H:%M:%S')
R = Rectangle(-6,36,30,46)
FLOAT_LIST=Float_Manager.FloatSelector(None, TI, R)
datetimelist = [f.time for f in FLOAT_LIST]
Coupled_List = TL.couple_with(datetimelist)
#  -----------------  init ----------------------

write_files_for_profiling('jobProfiler.sh', Coupled_List, FLOAT_LIST)
#execution of jobProfiler


# Coupled_List=[]
# for ir, req in enumerate(TL.getOwnList()):
#     LIST_of_REQ=[]
#     for f in FLOAT_LIST:
#         if (f.time >= req.starttime) & (f.time <= req.endtime) :
#             LIST_of_REQ.append(f)
#     if (len(LIST_of_REQ) >0 ): Coupled_List.append((TL.Timelist[ir],LIST_of_REQ))



