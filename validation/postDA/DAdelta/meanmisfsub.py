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
    Calculates mean misf for subbasins
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =True,
                                help = ''' Directory with misfit file 
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

TL_DIR = TimeList.fromfilenames(TI,RST_DIR,'*chl_mis.nc',  \
                dateformat='%Y%m%d',prefix='')

reset = False

masktypeLIST = ['open','coast','everywhere']
maskarea = {
    'open': TheMask.mask_at_level(200),
    'everywhere': TheMask.mask_at_level(0),
    'coast': TheMask.mask_at_level(0) & (TheMask.mask_at_level(200) == False),
}

misfTime = []
misfDate = []
varname = 'misfchl'
for itime,DAdate in enumerate(TL_DIR.Timelist):
    date8 = DAdate.strftime('%Y%m%d')
    print date8
    # filemisfit = TL_DIR.filelist[itime] + '/' \
    #             + date8 + '.chl_mis.nc'
    filemisfit = TL_DIR.filelist[itime]
    FM = NC.Dataset(filemisfit,'r')
    chlmis = np.array(FM.variables[varname])
    maskmis = (chlmis>1.e+19) | (chlmis<-900)
    mismean = {}
    for masktype in masktypeLIST:
        mismean[masktype] = np.zeros(Nsub)
        for isub,sub in enumerate(OGS.P):
            masksub = (maskmis==False) & (maskarea[masktype]==True) & (SUB[sub.name]==True)
            mismean[masktype][isub] = np.nanmean(chlmis[masksub])
    misfTime.append(mismean)
    misfDate.append(DAdate)


fileout = OUT_DIR + '/meanmisf_Time.npy'
np.save(fileout,misfTime)
fileout = OUT_DIR + '/meanmisf_Date.npy'
np.save(fileout,misfDate)











