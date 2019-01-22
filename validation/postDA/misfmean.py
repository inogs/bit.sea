import numpy as np
from commons.dataextractor import DataExtractor
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
from commons import timerequestors
from commons.season import season
import pickle

RUNname = 'RUN_TEST'
RUNname = 'RUN_2017Dard'
RUNname = 'TRANSITION'
RUNname = 'RUN_2017Dard_04'


DIRscratch = '/gpfs/scratch/userexternal/ateruzzi/'
DIRscratchgb = '/gpfs/scratch/userexternal/gbolzon0/'
DICTindir = {
        'TRANSITION' : DIRscratch + '/TRANS_24_for_comparison/DA_misf/',
        'RUN_TEST' : DIRscratch + '/TEST_2017aprjun/wrkdir/MODEL/DA__FREQ_1/',
        'RUN_2017Dard' : DIRscratchgb + 'OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/DA__FREQ_1/',
        'RUN_2017Dard_04' : DIRscratchgb + 'OPEN_BOUNDARY/TEST_04/wrkdir/MODEL/DA__FREQ_1/',
        }

DICToutdir = {
        'TRANSITION' : DIRscratch + '/ELAB_CFR_24/MISF_DA/TRANSITION/',
        'RUN_TEST' : DIRscratch + '/ELAB_CFR_24/MISF_DA/RUN_TEST/',
        'RUN_2017Dard' : DIRscratch + '/ELAB_CFR_24/MISF_DA/RUN_2017Dard/',
        'RUN_2017Dard_04' : DIRscratch + '/ELAB_CFR_24/MISF_DA/RUN_2017Dard_04/',
        }

INDIR = DICTindir[RUNname]
OUTdir = DICToutdir[RUNname]


maskfile = '/gpfs/scratch/userexternal/ateruzzi/TEST_2017aprjun/wrkdir/MODEL/meshmask.nc'
TheMask = Mask(maskfile)
_,jpj,jpi = TheMask.shape
mask200 = TheMask.mask_at_level(200)

Nsub = len(OGS.P.basin_list)
SUB = {}
SUB['med'] = np.zeros((jpj,jpi),np.bool)
for sub in OGS.Pred:
    print sub.name
    sbmask = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name] = sbmask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]


StartDate = '20170101'
EndDate   = '20180101'

TI = TimeInterval(StartDate,EndDate,'%Y%m%d')
misf_TL = TimeList.fromfilenames(TI, INDIR  ,"*chl_mis.nc", \
            prefix="", dateformat="%Y%m%d")

Ntime = misf_TL.nTimes

S=season()
misfmean = {}
misfmean['open'] = np.zeros((Nsub,Ntime))
misfmean['everywhere'] = np.zeros((Nsub,Ntime))
RMS = {}
RMS['open'] = np.zeros((Nsub,Ntime))
RMS['everywhere'] = np.zeros((Nsub,Ntime))

RMSseas = {}
RMSseas['win'] = np.zeros((2,Nsub))
RMSseas['spr'] = np.zeros((2,Nsub))
RMSseas['sum'] = np.zeros((2,Nsub))
RMSseas['aut'] = np.zeros((2,Nsub))


for ifile,misffile in enumerate(misf_TL.filelist):
    print misffile
    De = DataExtractor(TheMask,filename=misffile,varname='misfchl',dimvar=2)
    misf = De.values
    misf[misf>10e+10] = np.nan

    for isub,sub in enumerate(OGS.P):
        misfmean['everywhere'][isub,ifile] = np.nanmean(misf[SUB[sub.name]])
        RMS['everywhere'][isub,ifile] = (np.nanmean((misf[SUB[sub.name]])**2))**.5
        masksubopen = SUB[sub.name] & mask200
        misfmean['open'][isub,ifile] = np.nanmean(misf[masksubopen])
        RMS['open'][isub,ifile] = (np.nanmean((misf[masksubopen])**2))**.5

        iseas = 0
        seasReq = timerequestors.Season_req(2017,iseas,S)
        ii,_ = misf_TL.select(seasReq)
        RMSseas['win'][0,isub] = np.nanmean(RMS['everywhere'][isub,ii])
        RMSseas['win'][1,isub] = np.nanmean(RMS['open'][isub,ii])

        iseas = 1
        seasReq = timerequestors.Season_req(2017,iseas,S)
        ii,_ = misf_TL.select(seasReq)
        RMSseas['spr'][0,isub] = np.nanmean(RMS['everywhere'][isub,ii])
        RMSseas['spr'][1,isub] = np.nanmean(RMS['open'][isub,ii])

        iseas = 2
        seasReq = timerequestors.Season_req(2017,iseas,S)
        ii,_ = misf_TL.select(seasReq)
        RMSseas['sum'][0,isub] = np.nanmean(RMS['everywhere'][isub,ii])
        RMSseas['sum'][1,isub] = np.nanmean(RMS['open'][isub,ii])

        iseas = 3
        seasReq = timerequestors.Season_req(2017,iseas,S)
        ii,_ = misf_TL.select(seasReq)
        RMSseas['aut'][0,isub] = np.nanmean(RMS['everywhere'][isub,ii])
        RMSseas['aut'][1,isub] = np.nanmean(RMS['open'][isub,ii])



np.savetxt(OUTdir + '/subRMSDsum.txt',RMSseas['sum'])
np.savetxt(OUTdir + '/subRMSDspr.txt',RMSseas['spr'])
np.savetxt(OUTdir + '/subRMSDaut.txt',RMSseas['aut'])
np.savetxt(OUTdir + '/subRMSDwin.txt',RMSseas['win'])

LISTmisf = [0,0,0]

LISTmisf[0] = misf_TL.Timelist
LISTmisf[1] = misfmean
LISTmisf[2] = RMS


fileout = OUTdir + '/meanmisf.pkl'
fid = open(fileout,'wb')
pickle.dump(LISTmisf,fid)
fid.close()



