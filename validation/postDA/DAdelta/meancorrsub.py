import argparse
import numpy as np
import netCDF4 as NC

from commons.dataextractor import DataExtractor
# from commons.layer import Layer
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
# from layer_integral.mapbuilder import MapBuilder


def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates mean corr for subbasins
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =True,
                                help = ''' Directory with corr file 
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
Time__end="19990201"
Time__end="20500901"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")

TL_DIR = TimeList.fromfilenames(TI,RST_DIR,'*_corr.nc',  \
                dateformat='%Y%m%d',prefix='')

reset = False

masktypeLIST = ['open','coast','everywhere']
maskarea = {
    'open': TheMask.mask_at_level(200),
    'everywhere': TheMask.mask_at_level(0),
    'coast': TheMask.mask_at_level(0) & (TheMask.mask_at_level(200) == False),
}

corrTime = []
corrDate = []
varname = 'chl'
for itime,DAdate in enumerate(TL_DIR.Timelist):
    date8 = DAdate.strftime('%Y%m%d')
    print date8
    #filecorr = TL_DIR.filelist[itime] + '/' \
    #            + date8 + '_corr.nc'
    filecorr = TL_DIR.filelist[itime]
    FM = NC.Dataset(filecorr,'r')
    chlcorr = np.array(FM.variables[varname])
    maskcorr = (chlcorr>1.e+19) | (chlcorr==0)
    corrmean = {}
    for masktype in masktypeLIST:
        corrmean[masktype] = np.zeros(Nsub)
        for isub,sub in enumerate(OGS.P):
            masksub = (maskcorr==False) & (maskarea[masktype]==True) & (SUB[sub.name]==True)
            corrmean[masktype][isub] = np.nanmean(chlcorr[masksub])
    corrTime.append(corrmean)
    corrDate.append(DAdate)


fileout = OUT_DIR + '/meancorr_Time.npy'
np.save(fileout,corrTime)
fileout = OUT_DIR + '/meancorr_Date.npy'
np.save(fileout,corrDate)

