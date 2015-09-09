import numpy as np
import datetime

from region import Region, Rectangle

class Time_Interval():
    def __init__(self,starttime="19500101", endtime="21000101", dateformat='%Y%m%d'):
        self.starttime=datetime.datetime.strptime(starttime,dateformat)
        self.end__time=datetime.datetime.strptime(endtime  ,dateformat)

class Bio_Float():
    def __init__(self,lon,lat,time,filename):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.filename = filename

def FloatSelector(var, T, region):
    '''
    Arguments:
       var is a string indicating variable, 
          if var is None, no selection is done about variable
       T is as Time_Interval istance
       region is a region istance
    '''
    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )
    
    FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Index.txt"
    INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",")
    nFiles=INDEX_FILE.size    
    SELECTED=[]
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')
        Condition = var in available_params
        if var is None : Condition = True
             
        if Condition:
            if (float_time >= T.starttime) & (float_time <T.end__time):
                if region.is_inside(lon, lat):
                    SELECTED.append(Bio_Float(lon,lat,float_time,filename))
                      
    return SELECTED

if __name__ == '__main__':
    var = 'NITRATE'
    TI = Time_Interval('20150520','20150830','%Y%m%d')
    R = Rectangle(-6,36,30,46)
    
    FLOAT_LIST=FloatSelector(var, TI, R)

  
            
