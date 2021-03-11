import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Executed usually just after float download, creates Float_Index.0.txt file,
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = 'e.g. Med_floats.txt')

    parser.add_argument(   '--output_float_indexer','-o',
                                type = str,
                                required = True,
                                help = '''float indexer rough file as
                                /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_BIO/Float_Indexer.0.txt''')


    return parser.parse_args()

args = argument()
import numpy as np
from datetime import datetime

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


input_file="/gpfs/scratch/userexternal/gbolzon0/V7C/Med_floats.txt"
CORIOLIS_INDEX_FILE=np.loadtxt(args.inputfile,dtype=mydtype, delimiter=",")

nFiles = len(CORIOLIS_INDEX_FILE)


Float_Index_dtype= np.dtype([
          ('file_name','S200'),
          ('lat',np.float64),
          ('lon',np.float64),
          ('time','S17'),
          ('parameters','S200'),
          ('parameter_data_mode','S100')] )


INDEX_FILE = np.zeros((nFiles,), dtype=Float_Index_dtype)

INDEX_FILE['lon']                 = CORIOLIS_INDEX_FILE['longitude']
INDEX_FILE['lat']                 = CORIOLIS_INDEX_FILE['latitude']
INDEX_FILE['parameters']          = CORIOLIS_INDEX_FILE['parameters']
INDEX_FILE['parameter_data_mode'] = CORIOLIS_INDEX_FILE['parameter_data_mode']

for ifile in range(nFiles):
    filename_coriolis=CORIOLIS_INDEX_FILE[ifile]['file_name']
    _,wmo,_,basename = filename_coriolis.rsplit('/')
    d=datetime.strptime(CORIOLIS_INDEX_FILE[ifile]['date'],'%Y%m%d%H%M%S')
    INDEX_FILE[ifile]['file_name'] = wmo + "/" + basename
    INDEX_FILE[ifile]['time'] = d.strftime('%Y%m%d-%H:%M:%S')

np.savetxt(args.output_float_indexer, INDEX_FILE, fmt="%s,%f,%f,%s,%s,%s")
    