import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')
    parser.add_argument(   '--ingvmask', '-M',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')   
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')
    parser.add_argument(   '--prefix', '-p',
                                type = str,
                                default = "mfs_efs2_20171003_",
                                required = False,
                                help = ''' Typical prefix of ingv file, the important is the length''')
    parser.add_argument(   '--dateformat', '-f',
                                type = str,
                                default = "%Y%m%d",
                                required = False,
                                help = ''' dateformat''')
    parser.add_argument(   '--filepattern', '-l',
                                type = str,
                                default = "*nc",
                                required = False,
                                help = ''' file pattern''')    



    return parser.parse_args()

args = argument()
from commons.utils import addsep
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.Timelist import TimeInterval, TimeList
from commons import netcdf4


INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR = addsep(args.outdir)

MaskINGV   = Mask(args.ingvmask)
TheMask    = Mask(args.maskfile)

JPK, JPJ, JPI = MaskINGV.shape
jpk, jpj, jpi = TheMask.shape

lon_cut = JPI - jpi
 

TI = TimeInterval("1950","2050",'%Y')
TL = TimeList.fromfilenames(TI, INPUTDIR, args.filepattern, prefix=args.prefix, dateformat=args.dateformat)

for iFrame, ingv_file in enumerate(TL.filelist):
    outfile = OUTPUTDIR + TL.Timelist[iFrame].strftime("ave.%Y%m%d-12:00:00.votemper.nc")
    T = DataExtractor(MaskINGV,ingv_file,'votemper').values
    netcdf4.write_3d_file(T[:jpk,:,lon_cut:], 'votemper', outfile, TheMask)
    
    outfile = OUTPUTDIR + TL.Timelist[iFrame].strftime("ave.%Y%m%d-12:00:00.vosaline.nc")
    S = DataExtractor(MaskINGV,ingv_file,'vosaline').values
    netcdf4.write_3d_file(S[:jpk,:,lon_cut:], 'vosaline', outfile, TheMask)
    
    

