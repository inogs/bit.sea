import numpy as np
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
import pickle


RUNname = 'RUN_TEST'
RUNname = 'RUN_2017Dard'
RUNname = 'TRANSITION'


OUTdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/RST_DA/' + RUNname

DICTdir_b = {
    'TRANSITION': '/gpfs/scratch/userexternal/ateruzzi/' + \
                  'TRANS_24_for_comparison/DA_RST/',
    'RUN_TEST': '/gpfs/scratch/userexternal/ateruzzi/' + \
                'TEST_2017aprjun/wrkdir/MODEL/DA__FREQ_1/',
    'RUN_2017Dard': '/gpfs/scratch/userexternal/gbolzon0/' + \
                '/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/DA__FREQ_1/',
    }
DICTdir_a = {
    'TRANSITION': '/gpfs/scratch/userexternal/ateruzzi/TRANS_24_for_comparison/DA_RST/',
    'RUN_TEST': '/gpfs/scratch/userexternal/ateruzzi/' + \
                'TEST_2017aprjun/wrkdir/MODEL/RESTARTS/',
    'RUN_2017Dard': '/gpfs/scratch/userexternal/gbolzon0/' + \
                '/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/RESTARTS/',
    }
DICTfilename = {
    'TRANSITION': ['before','after.'],
    'RUN_TEST' : ['',''],
    'RUN_2017Dard' : ['',''],
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
EndDate   = '20180101'


varlist = ['P1l','P2l','P3l','P4l']

chlbef = {}
chlaft = {}
dates = []

for var in varlist:
    print '----------- ' + var + '-----------------' 
    TI = TimeInterval(StartDate,EndDate,'%Y%m%d')
    TL_b = TimeList.fromfilenames(TI,DICTdir_b[RUNname],
                    "RST*" + DICTfilename[RUNname][0] + "*00:*" + var  + "*nc",
                    prefix='RST.' + DICTfilename[RUNname][0])
    TL_a = TimeList.fromfilenames(TI,DICTdir_a[RUNname],
                    "RST*" + DICTfilename[RUNname][1] + "*00:*" + var  + "*nc",
                    prefix='RST.' + DICTfilename[RUNname][1])
    Ntime = TL_b.nTimes

    chlbef[var] = np.zeros((Nsub,Ntime,2))
    chlaft[var] = np.zeros((Nsub,Ntime,2))
    

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
        if var==varlist[0]: dates.append(dateRST)
        fileRSTa = TL_a.filelist[ind_a]
        print fileRSTa
        De = DataExtractor(TheMask,filename=fileRSTa,varname='TRN' + var,dimvar=3)
        chl_a = De.filled_values[0,:,:]

        for isub,sub in enumerate(OGS.P):
            chlbef[var][isub,ii,0] = np.nanmean(chl_b[SUB[sub.name]])
            chlbef[var][isub,ii,1] = np.nanmean(chl_b[SUB[sub.name] & mask200])
            chlaft[var][isub,ii,0] = np.nanmean(chl_a[SUB[sub.name]])
            chlaft[var][isub,ii,1] = np.nanmean(chl_a[SUB[sub.name] & mask200])

chltotbef =  {}
chltotaft =  {}

chltotbef['everywhere'] = np.zeros((Nsub,Ntime))
chltotbef['open'] = np.zeros((Nsub,Ntime))
chltotaft['everywhere'] = np.zeros((Nsub,Ntime))
chltotaft['open'] = np.zeros((Nsub,Ntime))

for var in varlist:
    chltotbef['everywhere'] = chltotbef['everywhere'] + chlbef[var][:,:,0]
    chltotbef['open'] = chltotbef['open'] + chlbef[var][:,:,1]
    chltotaft['everywhere'] = chltotaft['everywhere'] + chlaft[var][:,:,0]
    chltotaft['open'] = chltotaft['open'] + chlaft[var][:,:,1]


LISTb = [0,0]
LISTa = [0,0]

LISTb[0] = TL_b.Timelist
LISTb[1] = chltotbef
LISTa[0] = dates
LISTa[1] = chltotaft


fileout = OUTdir + '/meanRSTbefore.pkl'
fid = open(fileout,'wb')
pickle.dump(LISTb,fid)
fid.close()

fileout = OUTdir + '/meanRSTafter.pkl'
fid = open(fileout,'wb')
pickle.dump(LISTa,fid)
fid.close()


