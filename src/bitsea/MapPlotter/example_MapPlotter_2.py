#!/usr/bin/env python
#
# Example of how to use the MapPlotter class
# from already loaded data.
#
# Arnau Miro, OGS (2019)
from __future__ import print_function

import numpy as np
import MapPlotter as mp

from commons import netcdf4 as nc4
from commons.mask import Mask

## DATA PATHS ##
fname    = '/media/internal/disk2TiB/data/MEDICANE/AVE_PHYS/ave.20180925-12:00:00.votemper.nc'
varname  = 'votemper'
maskfile = '/media/internal/disk2TiB/data/MEDICANE/meshmask.nc'
idepth   = 0
outfile  = 'example_MapPlotter_2.png'
outdpi   = 300


# Load data
mask     = Mask(maskfile,zlevelsvar="nav_lev",ylevelsmatvar="gphit", xlevelsmatvar="glamt")
print(mask.xlevels.shape,mask.ylevels.shape)
data     = nc4.readfile(fname,varname)[0,:,:]
data[data > 1e10] = np.nan # Handle the mask
print(data.shape)

# Instance MapPlotter class
plotter = mp.MapPlotter(projection='PlateCarree')

# Create basic parameters dictionary
params  = plotter.defaultParams()
# DPI
params['dpi']      = 100
# Limits for the plot
params['xlim']     = [-6, 37]
params['ylim']     = [30, 46]
# Which features need to be plotted?
params['features'] = ['coastline','continents','rivers','image']
params['img']      = 'https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73726/world.topo.bathy.200406.3x5400x2700.png'
# A bit of formatting on the title and axis
params['title']    = ['Temperature Map',{'weight':'bold','style':'italic'}]
params['xlabel']   = ['Longitude (deg)']
params['ylabel']   = ['Latitude (deg)']
# A bit of formatting on the colorbar
params['cmap']     = 'jet'
params['bounds']   = [20,25]
params['label']    = {'label':'Temperature (deg C)','weight':None,'style':None}
print(params)

# Plot
plotter.plot(mask.xlevels,mask.ylevels,data,params=params)

# Save and show
plotter.save(outfile,dpi=outdpi)
plotter.show()