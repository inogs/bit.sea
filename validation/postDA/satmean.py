import numpy as np
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
# from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
import pickle
import Sat.SatManager as Sat


INDIR = '/gpfs/scratch/userexternal/ateruzzi/SATELLITE_TRANS/'


OUTdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/SAT/'


maskfile = '/gpfs/scratch/userexternal/ateruzzi/TEST_2017aprjun/wrkdir/MODEL/meshmask.nc'
TheMask = Mask(maskfile)
_,jpj,jpi = TheMask.shape
mask200 = TheMask.mask_at_level(200)

Nsub = len(OGS.P.basin_list)
SUB = {}
SUB['med'] = np.zeros((jpj,jpi),np.bool)
npointSub = {}
npointSub['everywhere'] = np.zeros(Nsub)
npointSub['open'] = np.zeros(Nsub)
for isub,sub in enumerate(OGS.Pred):
    print sub.name
    sbmask = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name] = sbmask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]
    npointSub['everywhere'][isub] = np.sum(SUB[sub.name]==True)
    npointSub['open'][isub] = np.sum((SUB[sub.name] & mask200)==True)

npointSub['everywhere'][-1] = np.sum(SUB['med']==True)
npointSub['open'][-1] = np.sum((SUB['med'] & mask200)==True)

fileout = OUTdir + '/pointsubbasin.pkl'
fid = open(fileout,'wb')
pickle.dump(npointSub,fid)
fid.close()

StartDate = '20170101'
EndDate   = '20180101'

TI = TimeInterval(StartDate,EndDate,'%Y%m%d')
sat_TL = TimeList.fromfilenames(TI, INDIR  ,"*.nc", \
            prefix="", dateformat="%Y%m%d")

Ntime = sat_TL.nTimes

satmean = {}
satmean['open'] = np.zeros((Nsub,Ntime))
satmean['everywhere'] = np.zeros((Nsub,Ntime))
satcoverage = {}
satcoverage['open'] = np.zeros((Nsub,Ntime))
satcoverage['everywhere'] = np.zeros((Nsub,Ntime))

for ifile,satfile in enumerate(sat_TL.filelist):
    print satfile
    chsat = Sat.readfromfile(satfile,var='CHL')
    chsat[chsat<0] = np.nan

    for isub,sub in enumerate(OGS.P):
        satmean['everywhere'][isub,ifile] = np.nanmean(chsat[SUB[sub.name]])
        satmean['open'][isub,ifile] = np.nanmean(chsat[SUB[sub.name] & mask200])
        satcoverage['everywhere'][isub,ifile] = np.sum(chsat[SUB[sub.name]]>0)
        satcoverage['open'][isub,ifile] = np.sum(chsat[SUB[sub.name] & mask200]>0)



LISTsat = [0,0,0]

LISTsat[0] = sat_TL.Timelist
LISTsat[1] = satmean
LISTsat[2] = satcoverage


fileout = OUTdir + '/meansat.pkl'
fid = open(fileout,'wb')
pickle.dump(LISTsat,fid)
fid.close()



