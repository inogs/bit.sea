import argparse
import numpy as np
import os

from commons.dataextractor import DataExtractor
from commons.layer import Layer
from commons.mask import Mask
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
from layer_integral.mapbuilder import MapBuilder


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

    parser.add_argument(   '--variable', '-v',
                                type = str,
                                required =True,
                                help = ''' Variable group (c, p, l, n)'''
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
vargroup = args.variable


TheMask = Mask(args.maskfile)
jpk,jpj,jpi = TheMask.shape

DA_layer = Layer(0,200)


Timestart="19990101"
Time__end="19991231"
Time__end="20500901"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")

TL_DIR = TimeList.fromfilenames(TI,RST_DIR,'*',  \
                dateformat='%Y%m%d',prefix='')


reset = False

for itime,DAdate in enumerate(TL_DIR.Timelist):
    date8 = DAdate.strftime('%Y%m%d')
    fileout = OUT_DIR + '/massbalances_' + vargroup + date8 + '.txt'
    if (os.path.exists(fileout)) and (not(reset)):
        print '  ... ' + date8 + ' already done'
        continue
    print date8
    deltaVar = np.zeros((jpj,jpi,4))
    IntBefore = np.zeros((jpj,jpi,4))
    IntAfter = np.zeros((jpj,jpi,4))
    totbal = np.zeros(3)
    totbal[:] = np.nan
    for iphyto in range(4):
        strphyto = np.str(iphyto + 1)
        varname = 'P' + strphyto + vargroup
        print ' ... ' + varname
        fileBefore = TL_DIR.filelist[itime] + '/RST.before.' \
                    + date8 + '-12:00:00.' \
                    + varname + '.nc'
        DeBEF = DataExtractor(TheMask,filename=fileBefore,varname=varname)
        IntBefore[:,:,iphyto] = MapBuilder.get_layer_integral(DeBEF,DA_layer)

        fileAfter = TL_DIR.filelist[itime] + '/RST.after.' \
                    + date8 + '-12:00:00.' \
                    + varname + '.nc'
        DeAFT = DataExtractor(TheMask,filename=fileAfter,varname='TRN' + varname)
        IntAfter[:,:,iphyto] = MapBuilder.get_layer_integral(DeAFT,DA_layer)
        
        deltaVar[:,:,iphyto] = IntAfter[:,:,iphyto] - IntBefore[:,:,iphyto]
    
    deltaVargroup = np.nansum(deltaVar,2)
    intBeforeVargroup = np.nansum(IntBefore,2)
    intAfterVargroup = np.nansum(IntAfter,2)
    deltaTOT = np.nansum(deltaVargroup*TheMask.area*10**-6)
    beforeTOT = np.nansum(intBeforeVargroup*TheMask.area*10**-6)
    afterTOT = np.nansum(intAfterVargroup*TheMask.area*10**-6)
    totbal[0] = deltaTOT
    totbal[1] = beforeTOT
    totbal[2] = afterTOT

    #np.savetxt(fileout,totbal)








