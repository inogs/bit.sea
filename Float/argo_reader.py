import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Read file argo, gives back file with mediterannean argo floats
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputfile',"-i",
                                type = str,
                                required = True,
                                help = 'input file argo txt, eg=')
    parser.add_argument(   '--outfile',"-o",
                                type = str,
                                required = True,
                                help = 'output file argo')
    return parser.parse_args()

args = argument()


import scipy.io.netcdf as NC
import numpy as np
import matplotlib.pyplot as pl

from basins.region import Rectangle
from basins import V2

input_file = args.inputfile
output_file  = args.outfile

mydtype= np.dtype([
          ('file_name','S200'),
	  ('date','S200'),
          ('latitude',np.float32),
          ('longitude',np.float32),
	  ('ocean','S10'),
	  ('profiler_type',np.int),
	  ('institution','S10'),
          ('parameters','S200'),
	  ('parameter_data_mode','S100'),
	  ('date_update','S200')] )

INDEX_FILE=np.loadtxt(input_file,dtype=mydtype, delimiter=",",ndmin=0,skiprows=9)


#file,date,latitude,longitude,ocean,profiler_type,institution,parameters,parameter_data_mode,date_update

#FIRST FILTER CORIOLIS

name_file=INDEX_FILE['file_name']
list_name=[] 
for name in name_file:           
    list_name.append(name.startswith('coriolis'))

Coriolis=INDEX_FILE[list_name]

#SECOND FILTER MEDITERRANEAN SEA

R = Rectangle(-6,36,30,46)
Mediterr=V2.med

lat_lon_list=[]
#for idx in Coriolis:  
#    ...:     lat_lon_list.append(R.is_inside(lat=idx['latitude'], lon=idx['longitude']))

for ele in Coriolis:
    lat_lon_list.append(Mediterr.is_inside(lat=ele['latitude'], lon=ele['longitude']))

MED_FLOAT=Coriolis[lat_lon_list]

np.savetxt(output_file, MED_FLOAT, fmt="%s,%s,%f,%f,%s,%d,%s,%s,%s,%s")


#THIRD FILTER BIOGEOCHEMICAL

#bioge_list=[]
#for bg in MED_FLOAT:
#    bioge_list.append(bg['parameters'] != 'PRES TEMP PSAL')

#BIOGE = MED_FLOAT[bioge_list]


#np.savetxt('bgc.txt', BIOGE, fmt="%s,%s,%f,%f,%s,%d,%s,%s,%s,%s")

