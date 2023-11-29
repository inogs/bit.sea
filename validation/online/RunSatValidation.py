import argparse

def argument():
    parser = argparse.ArgumentParser(description = 'run sat validation on N procs')
    parser.add_argument(   '--inputdir', '-i',
                                 type = str,
                                 required = True,
                                 help = 'Archive of P_l from model')

    parser.add_argument(   '--outputdir', '-o',
                                 type = str,
                                 required = True,
                                 help = 'Output dir of sat validation')


    parser.add_argument(   '--satdir', '-s',
                                 type = str,
                                 required = True,
                                 help = 'Daily satellite dir')

    parser.add_argument(   '--climdir', '-c',
                                 type = str,
                                 required = True,
                                 help = 'Clim dir')

    parser.add_argument(   '--maskfile', '-m',
                                 type = str,
                                 required = True,
                                 help = 'Maskfile')
    parser.add_argument(   '--var', '-v',
                                 type = str,
                                 required = True,
                                 choices = ['P_l','kd490'])
    
    parser.add_argument(   '--rundate', '-r',
                                 type = str,
                                 required = True,
                                 help = 'OPA_RUNDATE')

    return parser.parse_args()

args = argument()

from datetime import datetime,timedelta

from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
from commons.utils import addsep
import numpy as np
from SatValidation import SatValidation

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1



SAT_DAILY_DIR = addsep(args.satdir)
OUTDIR        = addsep(args.outputdir)
CLIM_DIR      = addsep(args.climdir)
ARCHIVEDIR    = addsep(args.inputdir)
OPA_RUNDATE   = args.rundate


TheMask=Mask(args.maskfile)

#basins
nSUB = len(OGS.P.basin_list)
jpk,jpj,jpi =TheMask.shape
mask200_2D = TheMask.mask_at_level(200.0)
dtype = [(sub.name, bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)

for sub in OGS.Pred:
    SUB[sub.name]  = SubMask(sub,maskobject=TheMask).mask_at_level(0)
    if 'atl' in sub.name: continue
    SUB['med'] = SUB['med'] | SUB[sub.name]

COASTNESS_LIST=['coast','open_sea','everywhere']

dtype = [(coast, bool) for coast in COASTNESS_LIST]
COASTNESS = np.ones((jpj,jpi),dtype=dtype)
COASTNESS['coast']     = ~mask200_2D
COASTNESS['open_sea']  =  mask200_2D

SUFFIX={'P_l'  : '_cmems_obs-oc_med_bgc-plankton_nrt_l3-multi-1km_P1D.nc',
        'kd490': '_cmems_obs-oc_med_bgc-transp_nrt_l3-multi-1km_P1D.nc'
        }

def climfilename(CLIM_DIR, month,var):
    if var=='P_l':
        return CLIM_DIR + '/ave.yyyy' + month + '15-00:00:00.nc'
    else:
        return None
        
delay=-15
 
opa_rundate=datetime.strptime(OPA_RUNDATE,'%Y%m%d')
for RUNDAY in range(1,8)[rank::nranks]:
    FC_RUNDATE = opa_rundate + timedelta(days=delay + RUNDAY)
    for fc_day in range(3):
        daydate = FC_RUNDATE + timedelta(days=fc_day)
        day = daydate.strftime('%Y%m%d')
        month = daydate.strftime('%m')
        avefile = "ave.%s-12:00:00.%s.nc" %(day, args.var)
        LOCAL_FC = "%s%s/%s" %(ARCHIVEDIR, FC_RUNDATE.strftime('%Y%m%d'),avefile)
        SAT_FILE = SAT_DAILY_DIR + day + SUFFIX[args.var]
        CLIM_FILE= climfilename(CLIM_DIR, month, args.var)
        f_name = OUTDIR + 'Validation_f' + str(fc_day+1) + '_' + FC_RUNDATE.strftime('%Y%m%d') + '_on_daily_Sat.' + day + '.nc'
        SatValidation(args.var, LOCAL_FC,SAT_FILE,CLIM_FILE,TheMask,f_name,SUB,COASTNESS_LIST,COASTNESS,nSUB)

    # PERSISTENCY  
    daydate = FC_RUNDATE + timedelta(days=-1)
    day = daydate.strftime('%Y%m%d')
    month = daydate.strftime('%m')
    avefile = "ave.%s-12:00:00.%s.nc" %(day, args.var)
    LOCAL_FC = "%s%s/%s" %(ARCHIVEDIR, FC_RUNDATE.strftime('%Y%m%d'),avefile)
    SATdate = FC_RUNDATE
    day_sat = SATdate.strftime('%Y%m%d')
    SAT_FILE=SAT_DAILY_DIR + day_sat + SUFFIX[args.var]
    CLIM_FILE= climfilename(CLIM_DIR, month, args.var)
    f_name = OUTDIR + 'Validation_pers' + '_' + FC_RUNDATE.strftime('%Y%m%d') + '_on_daily_Sat.' + day_sat + '.nc'
    SatValidation(args.var, LOCAL_FC,SAT_FILE,CLIM_FILE,TheMask,f_name,SUB,COASTNESS_LIST,COASTNESS,nSUB)

