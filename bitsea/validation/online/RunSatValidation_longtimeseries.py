import argparse

def argument():
     parser = argparse.ArgumentParser(description = 'run sat validation on HC - inputs ofr hystorical figure')
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
                                 help = 'Weekly satellite dir')

     parser.add_argument(   '--maskfile', '-m',
                                 type = str,
                                 required = True,
                                 help = 'Maskfile')

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



SAT_WEEKLY_DIR = addsep(args.satdir)
OUTDIR_HC = addsep(args.outputdir)
ARCHIVEDIR = addsep(args.inputdir)
OPA_RUNDATE   = args.rundate

TheMask=Mask(args.maskfile)

#basins
nSUB = len(OGS.P.basin_list)
_,jpj,jpi =TheMask.shape
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

SUFFIX = {
    'NRT': '_cmems_obs-oc_med_bgc-plankton_nrt_l3-multi-1km_P1D.nc',
    'DT' : '_cmems_obs-oc_med_bgc-plankton_my_l3-multi-1km_P1D.nc',

}

var = 'P_l'
delay = -7
CLIM_FILE = None

opa_rundate=datetime.strptime(OPA_RUNDATE,'%Y%m%d')
FC_RUNDATE = opa_rundate + timedelta(days=delay)
# FCdate = datetime.datetime.strptime(FC_RUNDATE,'%Y%m%d')

PREVIOUS_TUE_RUNDATE = (FC_RUNDATE - timedelta(days=7)).strftime('%Y%m%d')
HCdate = FC_RUNDATE - timedelta(days=8)
day = HCdate.strftime('%Y%m%d')
avefile = 'ave.' + day + '-12:00:00.P_l.nc'
LOCAL_HC = "%s%s/%s" %(ARCHIVEDIR, PREVIOUS_TUE_RUNDATE + '/POSTPROC/AVE_FREQ_1/ARCHIVE/',avefile)

SAT_WEEKLY_DIR_NRT = SAT_WEEKLY_DIR + '/NRT/WEEKLY_1_24/'
SATWEEKLY_FILE_NRT = SAT_WEEKLY_DIR_NRT + '/' + day + SUFFIX['NRT']

hc_name_NRT = OUTDIR_HC + '/Validation_hc_' + PREVIOUS_TUE_RUNDATE + '_on_weekly_SatNRT.' + day + '.nc'
SatValidation(var, LOCAL_HC,SATWEEKLY_FILE_NRT,CLIM_FILE,TheMask,hc_name_NRT,SUB,COASTNESS_LIST,COASTNESS,nSUB)

SAT_WEEKLY_DIR_DT = SAT_WEEKLY_DIR + '/DT/WEEKLY_1_24/'
SATWEEKLY_FILE_DT = SAT_WEEKLY_DIR_DT + '/' + day + SUFFIX['DT']

hc_name_DT = OUTDIR_HC + '/Validation_hc_' + PREVIOUS_TUE_RUNDATE + '_on_weekly_SatDT.' + day + '.nc'
SatValidation(var, LOCAL_HC,SATWEEKLY_FILE_NRT,CLIM_FILE,TheMask,hc_name_DT,SUB,COASTNESS_LIST,COASTNESS,nSUB)





