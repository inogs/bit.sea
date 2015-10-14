import numpy as np
import datetime
import pylab as pl
import os

from instruments.bio_float import BioFloat

class TimeInterval():
    def __init__(self, starttime="19500101", endtime="21000101", dateformat='%Y%m%d'):
        self.start_time = datetime.datetime.strptime(starttime, dateformat)
        self.end_time = datetime.datetime.strptime(endtime   ,dateformat)

    def contains(self, specific_time):
        if (specific_time >= self.start_time) and (specific_time < self.end_time):
            return True
        else:
            return False


def FloatSelector(var, T, region):
    '''
    Arguments:
       var is a string indicating variable, 
          if var is None, no selection is done about variable
       T is as TimeInterval istance
       region is a region istance
    '''
    mydtype= np.dtype([
              ('file_name','S200'),
              ('lat',np.float32),
              ('lon',np.float32),
              ('time','S17'),
              ('parameters','S200')] )
    
    FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Index.txt"
    #FloatIndexer="Float_Index.txt"
    INDEX_FILE=np.loadtxt(FloatIndexer,dtype=mydtype, delimiter=",",ndmin=1)
    nFiles=INDEX_FILE.size    
    selected = []
    for iFile in range(nFiles):
        timestr          = INDEX_FILE['time'][iFile]
        lon              = INDEX_FILE['lon' ][iFile]
        lat              = INDEX_FILE['lat' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d-%H:%M:%S')

        if var is None :
            VarCondition = True
        else:
            VarCondition = var in available_params

        if VarCondition:
            if T.contains(float_time) and region.is_inside(lon, lat):
                selected.append(BioFloat(lon,lat,float_time,filename,available_params))
                      
    return selected


