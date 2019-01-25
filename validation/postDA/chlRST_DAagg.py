import numpy as np
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
import pickle


RUNname = 'RUN_TEST'
RUNname = 'TRANSITION'
RUNname = 'RUN_2017Dard'


OUTdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/RST_DA/' + RUNname

DICTdir_b = {
    'TRANSITION': '/gpfs/scratch/userexternal/ateruzzi/' + \
                  'TRANS_24_for_comparison/DA_RSTaggbefore/',
    'RUN_2017Dard': '/gpfs/scratch/userexternal/ateruzzi/' + \
                  'RUN_2017Dard_for_comparison/DA_RSTaggbefore/',
    'RUN_TEST': '/gpfs/scratch/userexternal/ateruzzi/' + \
                'TEST_2017aprjun/wrkdir/POSTPROC/output/DA__FREQ_1/TMP/',
    }
DICTdir_a = {
    'TRANSITION': '/gpfs/scratch/userexternal/ateruzzi/' + \
                  'TRANS_24_for_comparison/DA_RSTaggafter/',
    'RUN_2017Dard': '/gpfs/scratch/userexternal/ateruzzi/' + \
                  'RUN_2017Dard_for_comparison/DA_RSTaggafter/',
    'RUN_TEST': '/gpfs/scratch/userexternal/ateruzzi/' + \
                'TEST_2017aprjun/wrkdir/POSTPROC/output/RESTARTS/TMP/',
    }
DICTfilename = {
    'TRANSITION': ['before','after.'],
    'RUN_2017Dard': ['',''],
    'RUN_TEST' : ['','']
    }

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
EndDate   = '20181201'


chlbef = {}
chlaft = {}
dates = []


var = 'P_l'
TI = TimeInterval(StartDate,EndDate,'%Y%m%d')
TL_b = TimeList.fromfilenames(TI,DICTdir_b[RUNname],
                "RST*" + DICTfilename[RUNname][0] + "*00:*" + var  + "*nc",
                prefix='RST.' + DICTfilename[RUNname][0])
TL_a = TimeList.fromfilenames(TI,DICTdir_a[RUNname],
                "RST*" + DICTfilename[RUNname][1] + "*00:*" + var  + "*nc",
                prefix='RST.' + DICTfilename[RUNname][1])
Ntime = TL_b.nTimes

chlbef['open'] = np.zeros((Nsub,Ntime))
chlaft['open'] = np.zeros((Nsub,Ntime))
chlbef['everywhere'] = np.zeros((Nsub,Ntime))
chlaft['everywhere'] = np.zeros((Nsub,Ntime))


for ii,fileRSTb in enumerate(TL_b.filelist):
    print fileRSTb
    De = DataExtractor(TheMask,filename=fileRSTb,varname=var,dimvar=3)
    chl_b = De.filled_values[0,:,:]
    dateRST = TL_b.Timelist[ii]
    ind_a = TL_a.find(dateRST)
    if not(dateRST==TL_a.Timelist[ind_a]):
        print 'Not found relevant "after" file'
        # continue
        break
    dates.append(dateRST)
    fileRSTa = TL_a.filelist[ind_a]
    print fileRSTa
    De = DataExtractor(TheMask,filename=fileRSTa,varname=var,dimvar=3)
    chl_a = De.filled_values[0,:,:]

    for isub,sub in enumerate(OGS.P):
        chlbef['everywhere'][isub,ii] = np.nanmean(chl_b[SUB[sub.name]])
        chlaft['everywhere'][isub,ii] = np.nanmean(chl_a[SUB[sub.name]])
        chlbef['open'][isub,ii] = np.nanmean(chl_b[SUB[sub.name] & mask200])
        chlaft['open'][isub,ii] = np.nanmean(chl_a[SUB[sub.name] & mask200])



LISTb = [0,0]
LISTa = [0,0]

LISTb[0] = TL_b.Timelist
LISTb[1] = chlbef
LISTa[0] = dates
LISTa[1] = chlaft


fileout = OUTdir + '/meanRSTbeforeagg.pkl'
fid = open(fileout,'wb')
pickle.dump(LISTb,fid)
fid.close()

fileout = OUTdir + '/meanRSTafteragg.pkl'
fid = open(fileout,'wb')
pickle.dump(LISTa,fid)
fid.close()


