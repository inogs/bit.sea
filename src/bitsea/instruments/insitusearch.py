import numpy as np
import os, datetime
import index_reader
dateformat17='%Y%m%d-%H:%M:%S'
# USER SELECTION ---------------------
DATESTART=datetime.datetime.strptime('20140101-00:00:00',dateformat17)
DATE__END=datetime.datetime.strptime('20150301-00:00:00',dateformat17)
var='DOX'
#-------------------------------------

mydtype= np.dtype([('catalog_id','S20'), 
          ('file_name','S200'),
          ('geospatial_lat_min',np.float32),
          ('geospatial_lat_max',np.float32),
          ('geospatial_lon_min',np.float32),
          ('geospatial_lon_max',np.float32),
          ('time_coverage_start','S19'),
          ('time_coverage_end','S19'),
          ('provider','S30'),
          ('date_update','S30'),
          ('data_mode','S1'),
          ('parameters','S200')] )


filename="/Users/gbolzon/Documents/OGS/COPERNICUS/medinsitu.hcmr.gr/index_monthly.txt"
A=np.loadtxt(filename,dtype=mydtype, delimiter=",", skiprows=6 )
nFiles=A.size





PARAMETERS=set()
for i in range(nFiles):
    mystr=str(A[i]['parameters']).strip()
    List=mystr.split(' ')
    if 'P' in List : print List
    for p in List :PARAMETERS.add(p)

VARLIST=['PHOS','SLCA','AMON','DOX1','DOX2','CPHL','NTRZ','NTRA','NTRI','PHPH']
#NTRZ:long_name = "nitrate + nitrite" ;



for i in range(nFiles):
    filename=os.path.basename(A['file_name'][i])
    parameters = A['parameters'][i]
    datestart=datetime.datetime.strptime(A[i]['time_coverage_start'],'%Y-%m-%dT%H:%M:%S')
    date_end=datetime.datetime.strptime(A[i]['time_coverage_end'  ],'%Y-%m-%dT%H:%M:%S')
    condition = ( var in parameters) & (filename[:2]=='MO') & (datestart > DATESTART) & (date_end < DATE__END) 
    if condition :
        print A[i]['file_name'], A[i]['parameters']

        