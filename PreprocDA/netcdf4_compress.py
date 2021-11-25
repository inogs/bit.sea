import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates NetCDF4 compressed files from NetCDF3.
    Works only for eofs files.
    Reduces data quality. 
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = '''Dir with NetCDF3 compressed files''')
    parser.add_argument(   '--maskfile','-m',
                                type = str,
                                required = True,
                                help = '''path of the mask file''')    
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                help = 'Dir with the NetCDF4 output compressed files ')
    return parser.parse_args()

args = argument()

from commons import netcdf4
import netCDF4 as NC
from commons.mask import Mask
import glob,os
from commons.utils import addsep

INPUTDIR = addsep(args.inputdir)
OUTDIR   = addsep(args.outdir)
TheMask = Mask(args.maskfile,loadtmask=False)
jk = TheMask.getDepthIndex(350.0)
filelist=glob.glob(INPUTDIR + "eof*nc")
filelist.sort()

for inputfile in filelist:
    outfile = OUTDIR + os.path.basename(inputfile)
    print(outfile)
    eva=netcdf4.readfile(inputfile, "eva")
    evc=netcdf4.readfile(inputfile, "evc")


    evc[:,jk:,:] = 0.0

    neof, nlev, nreg = evc.shape

    ncOUT = NC.Dataset(outfile,'w')
    ncOUT.createDimension("neof",neof)
    ncOUT.createDimension("nreg",nreg)
    ncOUT.createDimension("nlev",nlev)

    ncvar = ncOUT.createVariable("eva", 'f', ('neof','nreg'), zlib=True, fill_value=1.0e+20, complevel=9, least_significant_digit=2)
    ncvar[:]=eva
    ncvar = ncOUT.createVariable("evc", 'f', ('neof','nlev','nreg'), zlib=True, fill_value=1.0e+20,complevel=9, least_significant_digit=2)
    ncvar[:]=evc
    ncOUT.close()

