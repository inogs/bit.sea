# read_mask1x1_npy_do_netcdf.py

# read the file PYfile_Hlev_Hlevmask_mask1x1.npy created by readVAR_doMAP1DEG_13layer.py
# and save it as a netcdf file

import scipy.io.netcdf as NC
import glob
import os,sys
import numpy as np
curdir='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER'



[Hlev,Hlevmask,mask1x1]=np.load('PYfile_Hlev_Hlevmask_mask1x1.npy')


LON=np.arange(-8.5,35.501,1.0);  # equivalente del matlab -8.5:1:35.5
Xlon=np.arange(-8,35.001,1.0);
LAT=np.arange(29.5,46.501,1.0);  # equivalente del matlab 29.5:1:46.5
Ylat=np.arange(30,46.001,1.0);




#  scritttura del file NC\
nomefile = 'MASK1x1_13lev.nc'
outfile = nomefile
ncOUT=NC.netcdf_file(outfile,"w")
# write dimensions
ncOUT.createDimension('depth',13)
ncOUT.createDimension('lat',17)
ncOUT.createDimension('lon',44)

depth = ncOUT.createVariable('depth','f',('depth',))
depth[:] = [21.2,71.2,124.3, 177.2, 335.9, 721.5, 1210.7, 1758.9, 2293.9, 2798.4, 3294.5, 3752.3, 4269.3]
depth.units = 'depth of centre of the layers: 50m from 0 to 200, 200-500, 500m from 500 to 4500'

lon = ncOUT.createVariable('lon','f',('lon',))
lat = ncOUT.createVariable('lat','f',('lat',))
lon[:]=Xlon
lat[:]=Ylat
ncOUT.description = 'chain run with DA from 04-14 to 05-15 monthly: 1x1 degree , 13 layers'

#   save variable 
data = ncOUT.createVariable('mask1x1','f',('depth','lat','lon',))
data[:]=mask1x1
setattr(data,"missing_value",1.0e+20)
ncOUT.close()






