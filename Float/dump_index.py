import scipy.io.netcdf as NC
import datetime
import os,glob
# in the cronjob just after the download
#dump_index.py prints the index float file 
#e.g. lines as 
#/Users/gbolzon/Documents/OGS/COPERNICUS/bit.sea/MR6901765_024.nc 34.024883 24.519977 20150818-09:33:00 DOXY NITRATE CHLA PRES PSAL TEM


VARLIST=['DOXY','NITRATE','CHLA',  'PRES','PSAL','TEMP']


def file_header_content(filename,VARLIST):
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
    s="%s,%f,%f,%s," %(filename, lat, lon, Time.strftime('%Y%m%d-%H:%M:%S'))
    for var in VARLIST: 
        if var in ncIN.variables.keys():
            s=s+" " + var
    ncIN.close()
    return s

LOC="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/"
FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Index.txt"
DIRLIST=os.listdir(LOC)


LINES=[]
for DIR in DIRLIST:
    dirpath=LOC + DIR
    filenames = glob.glob(dirpath + "/*nc")
    for filename in filenames:
        if filename[-4:]!='D.nc':
            line=file_header_content(filename,VARLIST)
            if line is not None: LINES.append(line+"\n")

F = file(FloatIndexer,'w')
F.writelines(LINES)
F.close()

VARLIST=['DOXY','NO3','CHLA',  'PRES','PSAL','TEMP']
LOC="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_LOVBIO/"
FloatIndexer="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_LOVBIO/Float_Index.txt"
DIRLIST=os.listdir(LOC)


LINES=[]
for DIR in DIRLIST:
    dirpath=LOC + DIR
    filenames = glob.glob(dirpath + "/*nc")
    for filename in filenames:
        if filename[-4:]!='D.nc':
            line=file_header_content(filename,VARLIST)
            if line is not None: LINES.append(line+"\n")

F = file(FloatIndexer,'w')
F.writelines(LINES)
F.close()
