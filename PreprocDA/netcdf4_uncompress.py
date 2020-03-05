import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates NetCDF 3 files from NetCDF 4 compressed.
    Works only for eofs files.
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = '''NetCDF 4 compressed file''')
    parser.add_argument(   '--outfile','-o',
                                type = str,
                                required = True,
                                help = 'path of the NetCDF 3 output uncompressed file ')
    return parser.parse_args()

args = argument()

from commons import netcdf4
import scipy.io.netcdf as NC
eva=netcdf4.readfile(args.inputfile, "eva")
evc=netcdf4.readfile(args.inputfile, "evc")

neof, nlev, nreg = evc.shape

ncOUT = NC.netcdf_file(args.outfile,'w')
ncOUT.createDimension("neof",neof)
ncOUT.createDimension("nreg",nreg)
ncOUT.createDimension("nlev",nlev)

ncvar = ncOUT.createVariable("eva", 'f', ('neof','nreg'))
ncvar[:]=eva
ncvar = ncOUT.createVariable("evc", 'f', ('neof','nlev','nreg'))
ncvar[:]=evc
ncOUT.close()



