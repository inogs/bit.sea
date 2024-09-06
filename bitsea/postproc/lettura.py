import scipy.io.netcdf as NC
import numpy as np

import matplotlib.pyplot as pl

M=NC.netcdf_file("ave.20071005-12:00:00.N1p.nc","r")
N1p = M.variables['N1p'].data[0,:,:,:]
M.close()

A=N1p.copy()
A[A>1e+19]=np.NaN
pl.imshow(A[0,:,:])


pl.clim(0,.1)
pl.colorbar()

pl.gca().invert_yaxis()
pl.show()