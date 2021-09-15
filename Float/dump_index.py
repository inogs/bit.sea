import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Executed usually just after float download, creates Float_Index.0.txt file,
    the Database of floats, opening NetCDF files to take
    * lon,lat,time, parameters *
    If there is an operational Float_Indexer.txt of infos previously stored, dump_index.py
    takes them and opens only the just downloaded files.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/')
    parser.add_argument(   '--input_float_indexer','-f',
                                type = str,
                                required = False,
                                help = 'float indexer corrected file, like Float_Indexer.txt')
    parser.add_argument(   '--output_float_indexer','-o',
                                type = str,
                                required = True,
                                help = '''float indexer rough file as
                                /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Indexer.0.txt''')
    parser.add_argument(   '--type','-t',
                                type = str,
                                required = True,
                                choices = ['lov','coriolis','Float_opt', 'Float_opt_19', 'Float_opt_20','superfloat', 'static_superfloat'])

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
#6901765/MR6901765_024.nc 34.024883 24.519977 20150818-09:33:00 DOXY NITRATE CHLA PRES PSAL TEMP

NOW=datetime.datetime.now()
mydtype= np.dtype([
          ('file_name','S200'),
          ('lat',np.float64),
          ('lon',np.float64),
          ('time','S17'),
          ('parameters','S200')] )

FILELIST=[]
is_provided_indexer = False
if args.input_float_indexer is not None:
    if os.path.exists(args.input_float_indexer):
        INDEX_FILE=np.loadtxt(args.input_float_indexer,dtype=mydtype, delimiter=",",ndmin=1)
        FILELIST=INDEX_FILE['file_name'].tolist()
        is_provided_indexer = len(FILELIST) >0


if args.type=="coriolis":
    VARLIST=['DOXY','NITRATE','CHLA',  'PRES','PSAL','TEMP','PH_IN_SITU_TOTAL', 'BBP700','BBP532', 'DOWNWELLING_PAR','CDOM','DOWN_IRRADIANCE380'       ,'DOWN_IRRADIANCE412'       ,'DOWN_IRRADIANCE490' ]
if args.type=="lov":
    VARLIST=['DOXY','SR_NO3_ADJUSTED', 'CHLA',  'PRES','PSAL','TEMP','PH_IN_SITU_TOTAL', 'BBP700','BBP532', 'PAR'            ,'CDOM','DOWNWELLING_IRRADIANCE_380','DOWNWELLING_IRRADIANCE_412','DOWNWELLING_IRRADIANCE_490']
if args.type=='Float_opt':
    VARLIST=['PRES','PSAL','TEMP','PAR','CHLA', 'Ed_380','Ed_412','Ed_490']
if args.type=='Float_opt_19':
    VARLIST=['PRES', 'PAR','CHL','IRR_380','IRR_412','IRR_490']
if args.type=='Float_opt_20':
    VARLIST=['PRES', 'SALI','TEMP','BBP700']
if args.type=="superfloat":
    VARLIST=['DOXY','NITRATE','CHLA',  'PRES','PSAL','TEMP','PH_IN_SITU_TOTAL', 'BBP700','BBP532', 'DOWNWELLING_PAR','CDOM','DOWN_IRRADIANCE380'       ,'DOWN_IRRADIANCE412'       ,'DOWN_IRRADIANCE490' ]
if args.type=="static_superfloat":
    VARLIST=['DOXY','NITRATE','CHLA',  'PRES','PSAL','TEMP','PH_IN_SITU_TOTAL', 'BBP700','BBP532', 'PAR','CDOM','IRR_380' ,'IRR_412','IRR_490' ]


def file_header_content(filename,VARLIST, avail_params=None):
    '''
    it takes variable list
    Returns
    - a string like
        6901765/MR6901765_024.nc 34.024883 24.519977 20150818-09:33:00 DOXY NITRATE CHLA PRES PSAL TEMP
    - None in case of error
    '''
    try:
        ncIN = NC.netcdf_file(filename,'r')
    except:
        print("Not valid NetDCF file: " + filename)
        return

    lon=ncIN.variables['LONGITUDE'].data[0]
    lat=ncIN.variables['LATITUDE'].data[0]
    BadPosition = (lon > 90.) or (lon < -90.) or (lat > 90.) or (lat < -90.) 
    if BadPosition:
        print("Bad position in file : " + filename)
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
def get_sensor_list(wmo,LINES):
    for line in LINES:
        if wmo in line:
            d=StringIO(line)
            A=np.loadtxt(d,dtype=mydtype,delimiter=',')
            return str(A['parameters'])
    else:
        print(wmo + " not in CORIOLIS")
        return 'DOXY NITRATE CHLA PRES PSAL TEMP'


LOC=addsep(args.inputdir)
FloatIndexer=args.output_float_indexer
DIRLIST=os.listdir(LOC)
HERE=os.getcwd()
os.chdir(LOC)

LINES=[]
for DIR in DIRLIST:
    dirpath=DIR
    filenames = glob.glob(dirpath + "/*nc")
    filenames.sort()
    for filename in filenames:
        if filename[-4:]!='D.nc':
            if filename in FILELIST:
                ind=FILELIST.index(filename)
                timedist = NOW - datetime.datetime.strptime(INDEX_FILE['time'][ind][:8],"%Y%m%d")
                if timedist.days > 15:
                    line="%s,%f,%f,%s,%s" %(filename, INDEX_FILE['lat'][ind], INDEX_FILE['lon'][ind], INDEX_FILE['time'][ind], INDEX_FILE['parameters'][ind])
                    LINES.append(line+"\n")
                else:
                    line=file_header_content(filename,VARLIST,avail_params=None)
                    if line is not None:
                        if args.type=="lov": line = line.replace('SR_NO3_ADJUSTED','SR_NO3')
                        LINES.append(line+"\n")
            else:
                line=file_header_content(filename,VARLIST,avail_params=None)
                if line is not None:
                    if args.type=="lov": line = line.replace('SR_NO3_ADJUSTED','SR_NO3')
                    LINES.append(line+"\n")
                    if is_provided_indexer: print("added " + line)


F = file(FloatIndexer,'w')
F.writelines(LINES)
F.close()
os.chdir(HERE)

