import argparse
import numpy as np
import netCDF4 as NC
import os

from commons.dataextractor import DataExtractor
# from commons.layer import Layer
from commons.mask import Mask
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


Timestart="19990101"
Time__end="19990201"
Time__end="20500901"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")

TL_DIR = TimeList.fromfilenames(TI,RST_DIR,'*',  \
                dateformat='%Y%m%d',prefix='')

reset = False

misfTime = []
misfDate = []
varname = 'misfchl'
for itime,DAdate in enumerate(TL_DIR.Timelist):
    date8 = DAdate.strftime('%Y%m%d')
    print date8
    filemisfit = TL_DIR.filelist[itime] + '/' \
                + date8 + '.chl_mis.nc'
    FM = NC.Dataset(filemisfit,'r')
    chlmis = np.array(FM.variables[varname])
    maskmis = (chlmis>1.e+19) | (chlmis<-900)
    mismean = np.mean(chlmis[maskmis==False])
    misfTime.append(mismean)
    misfDate.append(DAdate)


fileout = OUT_DIR + '/meanmisf_Time.npy'
np.save(fileout,misfTime)
fileout = OUT_DIR + '/meanmisf_Date.npy'
np.save(fileout,misfDate)











