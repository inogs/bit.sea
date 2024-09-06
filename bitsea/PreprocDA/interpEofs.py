import numpy as np
import scipy.io.netcdf as NC
# import netCDF4
from commons.mask import Mask
import glob
import os
# from maskload import nav_lev as nav_lev16
from writeEOF import write_eof

maskfile24 = 'meshmask24.nc'
maskfile16 = 'meshmask16.nc'

Mask24 = Mask(maskfile24)
nav_lev24 = Mask24.zlevels
jpk24 = nav_lev24.shape[0]

Mask16 = Mask(maskfile16)
nav_lev16 = Mask16.zlevels
jpk16 = nav_lev16.shape[0]

varLIST = ['N3n']

DICTeoflev = {
    'P_l': 600,
    'N3n': 600,
}

DIREOF = 'EOF15sub_600/'

for var in varLIST:
    depEOF = DICTeoflev[var]
    ind16 = Mask16.getDepthIndex(depEOF)+1
    ind24 = Mask24.getDepthIndex(depEOF)+1 

    INDIR = DIREOF + '/' + var + '/'
    filelist = glob.glob(INDIR + '/eof.*.nc')
    for infile in filelist:
    # for infile in [filelist[0]]:
        print(infile)
        FF = NC.netcdf_file(infile,'r')
        eva = FF.variables['eva'].data.copy()
        evc16 = FF.variables['evc'].data.copy()
        neof, nlev, nreg = evc16.shape
        evc24 = np.zeros((neof,ind24+1,nreg))
        for ieof in range(neof):
            for ireg in range(nreg):
                evc24[ieof,:ind24+1,ireg] = np.interp(nav_lev24[:ind24+1], \
                                    nav_lev16[:ind16+1],evc16[ieof,:ind16+1,ireg], \
                                    left=evc16[ieof,0,ireg])
        fileout = var + '/' + os.path.basename(infile)
        print '... Writing ' + fileout
        write_eof(fileout,evc24,eva)

