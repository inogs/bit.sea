import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    interpolate variance field
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--indir', '-i',
                            type = str,
                            required = True,
                            help = 'Input variance'
                            )

    parser.add_argument(   '--maskfileorig', '-m',
                            type = str,
                            required = True,
                            help = 'Maskfile of the original resolution'
                            )

    parser.add_argument(   '--maskfilenew', '-n',
                            type = str,
                            required = True,
                            help = 'Maskfile of the final resolution'
                            )

    parser.add_argument(   '--prefix', '-p',
                            type = str,
                            required = True,
                            help = 'Maskfile of the final resolution'
                            )

    return parser.parse_args()

args = argument()

import numpy as np
import os
from commons.dataextractor import DataExtractor
from commons.interpolators import surf_interp_2d
from commons.mask import Mask
from commons import netcdf4
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.utils import addsep


INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)
PREFIX = args.prefix



MaskOrig = Mask(args.maskfileorig)
MaskNew = Mask(args.maskfilenew)


TL = TimeList.fromfilenames(None, INDIR, "*nc", prefix=PREFIX,  \
                            dateformat='%m')

masksea4 = MaskNew.mask[0,:,:]

for infile in TL.filelist:
    nomefile = os.path.basename(infile)
    outfile = OUTDIR + '/' + nomefile
    print outfile

    De = DataExtractor(MaskOrig,filename=infile,varname='variance',dimvar=2)
    varorig = De.filled_values

    # interpolation
    varnew = surf_interp_2d(MaskOrig,MaskNew,varorig)
    # if np.any(np.isnan(varnew[masksea4])):
    #     print np.sum(np.isnan(varnew[masksea4]))
    # varnew[~masksea4] = np.nan
    
    varnew[np.isnan(varnew)] = 1.e+20
    netcdf4.write_2d_file(varnew,'variance',outfile,MaskNew,fillValue=1.e+20)

