import numpy as np
import datetime

class Time_Interval():
    def __init__(self,starttime="19500101", endtime="21000101", dateformat='%Y%m%d'):
        self.starttime=datetime.datetime.strptime(starttime,dateformat)
        self.end__time=datetime.datetime.strptime(endtime  ,dateformat)

class Region():
    def __init__(self,lonmin,lonmax,latmin,latmax):
        self.lonmin = lonmin
        self.lonmax = lonmax
        self.latmin = latmin
        self.latmax = latmax

class Bio_Float():
    def __init__(self,lon,lat,time,filename):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename

def FloatSelector(var,T,Reg):
    '''
    Arguments:
       var is a string indicating variable
       T is as Time_Interval istance
       Reg is a Region istance
    '''
    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )
    
    INDEX_FILE=np.loadtxt("index.txt",dtype=mydtype, delimiter=",")
    nFiles=INDEX_FILE.size    
    SELECTED=[]
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')     
        if var in available_params:
            if (float_time >= T.starttime) & (float_time <T.end__time):
                if (lon >= Reg.lonmin) &  (lon <= Reg.lonmax) &  (lat >= Reg.latmin) &  (lat <= Reg.latmax):
                    SELECTED.append(Bio_Float(lon,lat,float_time,filename))
                      
    return SELECTED


var = 'NITRATE'
TI = Time_Interval('20150520','20150830','%Y%m%d')
R = Region(-6,36,30,46)

FLOAT_LIST=FloatSelector(var, TI, R)

  
            