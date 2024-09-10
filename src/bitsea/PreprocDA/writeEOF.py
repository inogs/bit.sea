import numpy as np
import scipy.io.netcdf as NC


def write_eof(outfile,eofp,eofa):
    neof, nlev, nreg = eofp.shape
    ncOUT   = NC.netcdf_file(outfile,"w")

    ncOUT.createDimension('nreg',nreg)
    ncOUT.createDimension('nlev',nlev)
    ncOUT.createDimension('neof',neof)

    ncvar = ncOUT.createVariable('eva','f',('neof','nreg'))
    ncvar[:] = eofa
    ncvar = ncOUT.createVariable('evc','f',('neof','nlev','nreg'))
    ncvar[:] = eofp

    ncOUT.close()


