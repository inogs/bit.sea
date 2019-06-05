import argparse
import numpy as np
import os

from commons.dataextractor import DataExtractor
from commons.layer import Layer
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
# from layer_integral.mapbuilder import MapBuilder


def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates mass delta du to DA
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =True,
                                help = ''' Directory with RST after and before 
                                (/marconi_scratch/userexternal/ateruzzi/ELAB_RA_COAST/DAdelta/DA_RST/)'''
                                )

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required =True,
                                help = ''' Maskfile meshmask.nc'''
                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required =True,
                                help = ''' output dir '''
                                )

    return parser.parse_args()


args = argument()

RST_DIR = addsep(args.inputdir)
OUT_DIR = addsep(args.outputdir)


TheMask = Mask(args.maskfile)
jpk,jpj,jpi = TheMask.shape

dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
Nsub = 0 
for sub in OGS.Pred:
    print sub
    sbmask         = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name]  = sbmask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]
    Nsub += 1
Nsub += 1


Timestart="19990101"
Time__end="19991231"
Time__end="20500901"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")

TL_DIR = TimeList.fromfilenames(TI,RST_DIR,'RST*P1l.nc',  \
                dateformat='%Y%m%d',prefix='RST.')


reset = False

masktypeLIST = ['open','coast','everywhere']
maskarea = {
    'open': TheMask.mask_at_level(200),
    'everywhere': TheMask.mask_at_level(0),
    'coast': TheMask.mask_at_level(0) & (TheMask.mask_at_level(200) == False),
}


deltaTime = []
deltaDate = []
for itime,DAdate in enumerate(TL_DIR.Timelist):
    date8 = DAdate.strftime('%Y%m%d')
#    if (os.path.exists(fileout)) and (not(reset)):
#        print '  ... ' + date8 + ' already done'
#        continue
    print date8
    deltaVar = np.zeros((jpj,jpi))
    chlBefore = np.zeros((jpj,jpi))
    chlAfter = np.zeros((jpj,jpi))
    # totbal = np.zeros(3)
    # totbal[:] = np.nan
    for iphyto in range(4):
        strphyto = np.str(iphyto + 1)
        varname = 'P' + strphyto + 'l'
        print ' ... ' + varname
      #  fileBefore = TL_DIR.filelist[itime] + '/RST.before.' \
      #              + date8 + '-12:00:00.' \
      #              + varname + '.nc'
      #  DeBEF = DataExtractor(TheMask,filename=fileBefore,varname=varname)
      #  chlBefore[:,:] = chlBefore[:,:] + DeBEF.filled_values[0,:,:]
        fileBefore = RST_DIR + '/RST.' \
                    + date8 + '-12:00:00.' \
                    + varname + '.nc'
        DeBEF = DataExtractor(TheMask,filename=fileBefore,varname=varname)
        chlBefore[:,:] = chlBefore[:,:] + DeBEF.filled_values[0,:,:]

      #  fileAfter = TL_DIR.filelist[itime] + '/RST.after.' \
      #              + date8 + '-12:00:00.' \
      #              + varname + '.nc'
        fileAfter = RST_DIR + '../RESTARTS/RST.' \
                    + date8 + '-12:00:00.' \
                    + varname + '.nc'
        DeAFT = DataExtractor(TheMask,filename=fileAfter,varname='TRN' + varname)
        chlAfter[:,:] = chlAfter[:,:] + DeAFT.filled_values[0,:,:]
    
    
    deltaVar[:,:] = chlAfter[:,:] - chlBefore[:,:]
    deltamean = {}
    for masktype in masktypeLIST:
        deltamean[masktype] = np.zeros(Nsub)
        for isub,sub in enumerate(OGS.P):
            masksub = (maskarea[masktype]==True) & (SUB[sub.name]==True)
            deltamean[masktype][isub] = np.nanmean(deltaVar[masksub])
    deltaTime.append(deltamean)
    deltaDate.append(DAdate)


fileout = OUT_DIR + '/meandelta_Time.npy'
np.save(fileout,deltaTime)
fileout = OUT_DIR + '/meandelta_Date.npy'
np.save(fileout,deltaDate)

    








