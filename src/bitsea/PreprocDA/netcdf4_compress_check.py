import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates NetCDF4 compressed files from NetCDF3.
    Works only for eofs files.
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = '''NetCDF 3 compressed file''')
    parser.add_argument(   '--maskfile','-m',
                                type = str,
                                required = True,
                                help = '''path of the mask file''')      
    parser.add_argument(   '--outfile','-o',
                                type = str,
                                required = True,
                                help = 'path of the NetCDF 4 output compressed file ')
    return parser.parse_args()

args = argument()

from bitsea.commons import netcdf4
import pylab as pl
from bitsea.commons.mask import Mask

TheMask=Mask(args.maskfile)

evc_old=netcdf4.readfile(args.inputfile,"evc")
evc_new=netcdf4.readfile(args.outfile,  "evc")

fig,ax=pl.subplots()

for ii in range(3):
    ax.plot(TheMask.zlevels[:60],evc_old[ii,:,33].T,'r',linewidth=3-ii)
    ax.plot(TheMask.zlevels[:60],evc_new[ii,:,33].T,'g',linewidth=3-ii)


fig.show()