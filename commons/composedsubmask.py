# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import os
import numpy as np
import netCDF4
from commons.mask import *
from basins.basin import ComposedBasin

class ComposedSubMask(Mask):
    """ Defines a submask starting from a ComposedBasin object and a NetCDF file or a Mask object
    """
    def __init__(self, basin, filename=None, maskobject=None, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t"):
        """ ComposedSubMask constructor

        Args:
            - *basin*: a ComposedBasin object.
            - *filename*: the path to a NetCDF file.
            - *maskobject*: a Mask object.

        Caveats:
            - *filename* and *maskobject* are mutually exclusive.
        """
        # If filename and maskobject are both defined or undefined
        if not ((filename is None) ^ (maskobject is None)):
            raise ValueError("You have to provide either a file name or a Mask object")
        elif not (maskobject is None):
            if not isinstance(maskobject, Mask):
                raise ValueError("maskobject must be an instance of Mask")
            base_mask = maskobject.mask
            self._mask = None
            self._shape = self._mask.shape
            self._xlevels = maskobject.xlevels
            self._ylevels = maskobject.ylevels
            self._zlevels = maskobject.zlevels
            self._dz = maskobject.dz
        elif not (filename is None):
            super(ComposedSubMask, self).__init__(filename, maskvarname, zlevelsvar, ylevelsmatvar, xlevelsmatvar, dzvarname)
            base_mask = self._mask
            self._mask = None
        if isinstance(basin, ComposedBasin):
            self._basin = basin
        else:
            raise ValueError("basin must be an instance of ComposedBasin")
        # Initialize composed mask
        dtype = [(b.name, np.bool) for b in self._basin]
        self._mask = np.zeros_like(base_mask, dtype=dtype)
        # For each basin
        for b in basin:
            # Build the prism mask
            prism = np.zeros_like(base_mask)
            for x in range(prism.shape[2]):
                for y in range(prism.shape[1]):
                    lon,lat = self.convert_x_y_to_lon_lat(x,y)
                    prism[:,y,x] = b.is_inside(lon,lat)
            # submask = original mask * prism mask
            self._mask[b.name] = base_mask * prism

    def save_as_netcdf(self, filename, maskvarname_suffix='_mask'):
        """Saves the submask data to a NetCDF file.

        Args:
            - *filename*: the path to the NetCDF file to save.
            - *descr*: a brief description of the data (optional, defaults to
              repr(self._basin) + " mask").
            - *maskvarname_suffix*: the suffix for the name of the variable to
              use inside the NetCDF file (optional, defaults to '_mask').
        """
        if os.path.exists(filename):
            # Open the file for appending
            netCDF_out = netCDF4.Dataset(filename, "a", format="NETCDF4")
        else:
            # Open the file for writing
            netCDF_out = netCDF4.Dataset(filename, "w", format="NETCDF4")

            # Add the spatial dimensions
            netCDF_out.createDimension('z', self._shape[0])
            netCDF_out.createDimension('y', self._shape[1])
            netCDF_out.createDimension('x', self._shape[2])

            # Add nav_lev data
            nav_lev = netCDF_out.createVariable('nav_lev', 'f4', ('z',))
            nav_lev[:] = self._zlevels

            # Create the nav_lat NetCDF variable
            nav_lat = netCDF_out.createVariable('nav_lat', 'f4', ('y', 'x'))
            nav_lat[:,:] = self._ylevels[:,:]

            # Create the nav_lon NetCDF variable
            nav_lon = netCDF_out.createVariable('nav_lon', 'f4', ('y', 'x'))
            nav_lon[:,:] = self._xlevels[:,:]

        # For each sub-basin
        for b in self._basin:
            # Set mask variable name
            maskvarname = b.name + maskvarname_suffix

            # Set description
            descr = repr(b) + " mask"

            # Create a variable to hold the data
            mask = netCDF_out.createVariable(maskvarname, 'u1', ('z', 'y', 'x'))

            # Assign the data to the NetCDF variable
            mask[:,:,:] = self._mask[b.name][:,:,:]

            # Write the description
            mask.descr = descr

        # Close the file
        netCDF_out.close()
