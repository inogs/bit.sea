import scipy.io.netcdf as NC
import datetime
import glob
# in the cronjob just after the download
#dump_index.py prints the index float file 
#e.g. lines as 
#/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/MR6901765_024.nc 34.024883 24.519977 20150818-09:33:00 DOXY NITRATE CHLA PRES PSAL TEM
# To be provided: input dir LOC and file for stdout


LOC="/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/"

filenames = glob.glob(LOC + "*nc")

VARLIST=['DOXY','NITRATE','CHLA',  'PRES','PSAL','TEMP']


def file_header_content(filename):
    ncIN = NC.netcdf_file(filename,'r')
    lon=ncIN.variables['LONGITUDE'].data[0]
    lat=ncIN.variables['LATITUDE'].data[0]
    
    ref  = ncIN.variables['REFERENCE_DATE_TIME'].data.tostring()
    juld = ncIN.variables['JULD'].data[0]
    d=datetime.datetime.strptime(ref,'%Y%m%d%H%M%S')
    Time =  d+datetime.timedelta(days=juld)
    s="%s,%f,%f,%s," %(filename, lat, lon, Time.strftime('%Y%m%d-%H:%M:%S'))
    for var in VARLIST: 
        if var in ncIN.variables.keys():
            s=s+" " + var
    ncIN.close()
    return s


for filename in filenames:
    print file_header_content(filename)