import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Reads the temporary json file obtainedn by maos_query.sh
    Dumps the wmo.txt file used by float harvesters.
    
    The expected practice is:
    
    ./maos_query.sh > $HOME_CHAIN/var/tmp/json_file
    if  [ $? -eq 0 ] ; then
        export WMO_FILE=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/FLOAT_LOVBIO/wmo.txt
        python json2table -j  $HOME_CHAIN/var/tmp/json_file -o $WMO_FILE 
    fi
    
    
    '''
    ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--jsonfile','-j',
                                type = str,
                                required = True,
                                help = 'path of the json file')
    parser.add_argument(   '--wmofile','-o',
                                type = str,
                                required = True,
                                help = 'path of the $DRES wmo file')

    return parser.parse_args()

args = argument()

import numpy as np
import sys, json
fid=open(args.jsonfile,'r') 
A=json.load(fid)  # list of dicts
fid.close()

nFloats= len(A)
FLOAT_dtype=[('id',np.int),('wmo','S20'),('id_type',np.int),('type','S20'),('nome_lov','S20'),('status','S1')]

TABLE=np.ones((nFloats,),dtype=FLOAT_dtype)


header="id_float |   wmo   | id_type |     type      |  nome_fs   | status\n" + "----------+---------+---------+---------------+------------+--------"
fmt="%10.f | %s | %d | %s | %s | %s"
for iFloat in range(nFloats):
    D = A[iFloat]
    
    TABLE[iFloat]['id']      =int(D['id_float'])
    TABLE[iFloat]['wmo']     =str(D['wmo'])
    TABLE[iFloat]['id_type'] =int(D['id_type'])
    TABLE[iFloat]['type']    =str(D['type'])
    TABLE[iFloat]['nome_lov'] =str(D['nome_lov'])
    TABLE[iFloat]['status']=  str(D['status'])
      

np.savetxt(args.wmofile, TABLE, fmt, header=header)
