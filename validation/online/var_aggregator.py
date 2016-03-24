import numpy as np
import scipy.io.netcdf as NC
import glob
import os
import argparse
from commons.utils import addsep

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''')
    
    
    parser.add_argument(   '--biodir',
                                type = str,
                                required = True,
                                help = '/some/path/MODEL/AVE_FREQ_1/')
    parser.add_argument(   '--physdir',
                                type = str,
                                required = True,
                                help = '/some/path/MODEL/AVE_FREQ_1/')

    parser.add_argument(   '--tmpdir', '-t',
                                type = str,
                                default = None,
                                required = True,
                                help = """ /some/path/POSTPROC/output/AVE_FREQ_1/TMP/ .  
                                Path to put files with aggregated variables for aveScan.py. 
                                """)
  
    return parser.parse_args()

def WriteTMPave(biofile,physfile, outfile):
            
    nc=NC.netcdf_file(biofile,"r");
    DIMS=nc.dimensions;    
    jpk = DIMS['depth']
    jpj = DIMS['lat'  ]
    jpi = DIMS['lon'  ]

    
    ncOUT=NC.netcdf_file(outfile,"w")
    setattr(ncOUT,"Convenctions","COARDS")
    setattr(ncOUT,"DateStart",nc.DateStart)
    setattr(ncOUT,"Date__End",nc.Date__End)
        
    
    for dimName,dimValue in DIMS.items():
        ncOUT.createDimension(dimName,dimValue)
    
    for var in ['lon','lat','depth']:
        ncvar=ncOUT.createVariable(var,'f',(var,))
        ncvar[:]=nc.variables[var].data
        setattr(ncvar,"actual_range",nc.variables[var].actual_range)
        setattr(ncvar,"units"       ,nc.variables[var].units)
    nc.close()        
    
    setattr(ncOUT.variables['lon'],"long_name","Longitude")    
    setattr(ncOUT.variables['lat'],"long_name","Latitude")

    
    for var in ['N1p','N3n','O2o']  :
        ncIN = NC.netcdf_file(biofile,"r")      
        ncvar=ncOUT.createVariable(var,'f',('time','depth','lat','lon'))
        ncvar[:]=ncIN.variables[var].data.copy()
        setattr(ncvar,"long_name",var)
        setattr(ncvar,"missing_value",1.e+20)
        ncIN.close()
    for var in ['votemper','vosaline']  :
        ncIN = NC.netcdf_file(physfile,"r")      
        ncvar=ncOUT.createVariable(var,'f',('time','depth','lat','lon'))
        ncvar[:]=ncIN.variables[var].data.copy()
        setattr(ncvar,"long_name",var)
        setattr(ncvar,"missing_value",1.e+20)
        ncIN.close()
    
    AGGREGATE_DICT={'P_i':['P1i','P2i','P3i','P4i']}
    for var in AGGREGATE_DICT.keys():
              
        ncvar=ncOUT.createVariable(var,'f',('time','depth','lat','lon'))
        junk = np.zeros((1,jpk,jpj,jpi),np.float32)
        for lvar in AGGREGATE_DICT[var]:
            ncIN = NC.netcdf_file(biofile,"r")     
            junk +=ncIN.variables[lvar].data.copy()
            ncIN.close()
        tmask= junk > 1.e+19
        junk[tmask] = 1.e+20
        ncvar[:]=junk    
        setattr(ncvar,"long_name",var)
        setattr(ncvar,"missing_value",1.e+20)
    
    ncOUT.close()    



try :
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks = comm.size 
except:
    rank   = 0
    nranks = 1

args = argument()
 

BIOAVEDIR     = addsep(args.biodir)
PHYS_DIR      = addsep(args.physdir)
PATH_NAME = BIOAVEDIR + "ave.*nc"
if rank==0 : 
    print "BIO_INPUT_DIR =", BIOAVEDIR


TMPOUTdir  = addsep(args.tmpdir)
if rank==0 : print "TMPOUTDIR= ", TMPOUTdir
os.system("mkdir -p " + TMPOUTdir)

archived_filelist=glob.glob(PATH_NAME)
archived_filelist.sort()

for filename in archived_filelist[rank::nranks]:
    dailyAve  = os.path.basename(filename)
    print "writing ", dailyAve
    daily_phys = PHYS_DIR + "T" + dailyAve[4:]
    outfile   = TMPOUTdir + dailyAve
    
        
    WriteTMPave(filename, daily_phys, outfile)
    

