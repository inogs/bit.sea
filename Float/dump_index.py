import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates Float_Index.txt files.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--coriolis',
                                type = str,
                                required = False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/",
                                help = 'directory of coriolis dataset')
    parser.add_argument(   '--lov',
                                type = str,
                                required = False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_LOVBIO/",
                                help = 'directory of lov dataset ')

    return parser.parse_args()

args = argument()


import scipy.io.netcdf as NC
import datetime
import os,glob
import numpy as np
from StringIO import StringIO
from commons.utils import addsep
# in the cronjob just after the download
#dump_index.py prints the index float file 
#e.g. lines as 
#/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/MR6901765_024.nc 34.024883 24.519977 20150818-09:33:00 DOXY NITRATE CHLA PRES PSAL TEM


VARLIST_COR=['DOXY','NITRATE','CHLA',  'PRES','PSAL','TEMP','PH_IN_SITU_TOTAL', 'BBP700','BBP532', 'DOWNWELLING_PAR','CDOM','DOWN_IRRADIANCE380'       ,'DOWN_IRRADIANCE412'       ,'DOWN_IRRADIANCE490' ]
VARLIST_LOV=['DOXY','SR_NO3', 'CHLA',  'PRES','PSAL','TEMP','PH_IN_SITU_TOTAL', 'BBP700','BBP532', 'PAR'            ,'CDOM','DOWNWELLING_IRRADIANCE_380','DOWNWELLING_IRRADIANCE_412','DOWNWELLING_IRRADIANCE_490']
VARLIST_OPT=['PRES','PSAL','TEMP','PAR','CHLA', 'Ed_380','Ed_412','Ed_490']

def file_header_content(filename,VARLIST, avail_params=None):
    '''
    it takes variable list
    '''
    ncIN = NC.netcdf_file(filename,'r')
    lon=ncIN.variables['LONGITUDE'].data[0]
    lat=ncIN.variables['LATITUDE'].data[0]
    BadPosition = (lon > 90.) or (lon < -90.) or (lat > 90.) or (lat < -90.) 
    if BadPosition:
        ncIN.close()
        return

    ref  = ncIN.variables['REFERENCE_DATE_TIME'].data.tostring()
    juld = ncIN.variables['JULD'].data[0]
    d=datetime.datetime.strptime(ref,'%Y%m%d%H%M%S')
    Time =  d+datetime.timedelta(days=juld)
    split_path=filename.rsplit(os.sep)
    wmo = split_path[-2]
    basename=split_path[-1]
    relative_name=wmo + "/" + basename
    s="%s,%f,%f,%s," %(relative_name, lat, lon, Time.strftime('%Y%m%d-%H:%M:%S'))
    if avail_params is None:
        for var in VARLIST: 
            if var in ncIN.variables.keys():
                s=s+" " + var
        ncIN.close()
    else:
        s=s+" " + avail_params
    return s

if True:
    LOC=addsep(args.coriolis) #"/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/"
    FloatIndexer=LOC + "Float_Index.0.txt"
    DIRLIST=os.listdir(LOC)


    LINES=[]
    for DIR in DIRLIST:
        dirpath=LOC + DIR
        filenames = glob.glob(dirpath + "/*nc")
        filenames.sort()
        for filename in filenames:
            if filename[-4:]!='D.nc':
                line=file_header_content(filename,VARLIST_COR,avail_params=None)
                if line is not None: LINES.append(line+"\n")

    F = file(FloatIndexer,'w')
    F.writelines(LINES)
    F.close()

    CORIOLIS_LINES=LINES[:]

if True:
    LOC=addsep(args.lov)
    FloatIndexer=LOC + "Float_Index.0.txt"
    DIRLIST=os.listdir(LOC)


    def get_sensor_list(wmo,LINES):
        mydtype= np.dtype([
                  ('file_name','S200'),
                  ('lat',np.float32),
                  ('lon',np.float32),
                  ('time','S17'),
                  ('parameters','S200')] )
        for line in LINES:
            if wmo in line:
                d=StringIO(line)
                A=np.loadtxt(d,dtype=mydtype,delimiter=',')
                return str(A['parameters'])
        else:
            print wmo + " not in CORIOLIS"
            return 'DOXY NITRATE CHLA PRES PSAL TEMP'



    LINES=[]
    for DIR in DIRLIST:
        wmo = DIR
        dirpath=LOC + DIR
        #sensors = get_sensor_list(wmo,CORIOLIS_LINES)
        filenames = glob.glob(dirpath + "/*nc")
        filenames.sort()
        for filename in filenames:
            line=file_header_content(filename,VARLIST_LOV,avail_params=None) #sensors.replace('NITRATE','SR_NO3'))
            if line is not None: LINES.append(line+"\n")

    F = file(FloatIndexer,'w')
    F.writelines(LINES)
    F.close()
