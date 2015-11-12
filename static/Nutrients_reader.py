import numpy as np
import scipy.io.netcdf as NC


filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Nutrients/Dataset_Med_Nutrients.nc"

ncIN=NC.netcdf_file(filename,'r')
DATA     = ncIN.variables['DATA'].data.copy()
UNITS    = ncIN.variables['UNITS'].data.copy()
VARIABLES= ncIN.variables['VARIABLES'].data.copy()
CRUISES  = ncIN.variables['Cruises'].data.copy()
ncIN.close()