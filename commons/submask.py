# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import os
import numpy as np
import netCDF4
from commons.mask import *
from commons.helpers import is_number
from basins.basin import Basin
from basins.basin import SimpleBasin
from basins.region import Rectangle

class SubMask(Mask):
    """ Defines a submask starting from a Basin object and a NetCDF file or a Mask object
    """
    def __init__(self, basin, filename=None, maskobject=None, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t"):
        """ SubMask constructor

        Args:
            - *basin*: a Basin object.
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
            self._mask = maskobject.mask
            self._shape = self._mask.shape
            self._xlevels = maskobject.xlevels
            self._ylevels = maskobject.ylevels
            self._zlevels = maskobject.zlevels
            self._dz = maskobject.dz
        elif not (filename is None):
            super(SubMask, self).__init__(filename, maskvarname, zlevelsvar, ylevelsmatvar, xlevelsmatvar, dzvarname)
        if isinstance(basin, Basin):
            self._basin = basin
        else:
            raise ValueError("basin must be an instance of Basin")
        # Build the prism mask
        prism = np.zeros_like(self._mask)
        for x in range(prism.shape[2]):
            for y in range(prism.shape[1]):
                lon,lat = self.convert_x_y_to_lon_lat(x,y)
                prism[:,y,x] = basin.is_inside(lon,lat)
        # submask = original mask * prism mask
        self._mask = self._mask * prism

    def save_as_netcdf(self, filename, descr=None, maskvarname='mask'):
        """Saves the submask data to a NetCDF file.

        Args:
            - *filename*: the path to the NetCDF file to save.
            - *descr*: a brief description of the data (optional, defaults to
              repr(self._basin) + " mask").
            - *maskvarname*: the name of the variable to use inside the NetCDF
              file (optional, defaults to 'mask').
        """
        # Set description
        if descr is None:
            descr = repr(self._basin) + " mask"

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

            # Prepare a zero-filled matrix capable of holding the nav_lat and nav_lon values
            zero_mat = np.zeros((self._shape[1], self._shape[2]), np.float32)

            # Create the nav_lat NetCDF variable
            nav_lat = netCDF_out.createVariable('nav_lat', 'f4', ('y', 'x'))
            nav_lat[:,:] = self._ylevels[:,np.newaxis] + zero_mat

            # Create the nav_lon NetCDF variable
            nav_lon = netCDF_out.createVariable('nav_lon', 'f4', ('y', 'x'))
            nav_lon[:,:] = self._xlevels[np.newaxis,:] + zero_mat

        # Create a variable to hold the data
        mask = netCDF_out.createVariable(maskvarname, 'u1', ('z', 'y', 'x'))

        # Assign the data to the NetCDF variable
        mask[:,:,:] = self._mask[:,:,:]

        # Write the description
        mask.descr = descr

        # Close the file
        netCDF_out.close()

    @staticmethod
    def from_square_cutting(mask, degrees, start_lon, start_lat):
        """
        Creates a list of SubMask objects from cutting a Mask object into
        square sections.  The cutting starts from a point and proceeds to cut
        along longitude then along latitude.

        ++++++++    ++++++++    s->+++++
        ++++++++ -> s->+++++ -> ssssssss
        s->+++++    ssssssss    ssssssss

        Args:
            - *mask*: a Mask object.
            - *degrees*: number of degrees for a sigle section side.
            - *start_lon*: starting longitude point.
            - *start_lat*: starting latitude point.

        Returns: a list of SubMasks mapping a single section each.
        """
        # Input validation
        if not isinstance(mask, Mask):
            raise ValueError("mask must be an instance of Mask")
        if not is_number(degrees):
            raise ValueError("degrees must be a number")
        if not is_number(start_lon):
            raise ValueError("start_lon must be a number")
        if not is_number(start_lat):
            raise ValueError("start_lat must be a number")
        # Get mask dimensions
        min_lon = min(mask.xlevels)
        max_lon = max(mask.xlevels)
        min_lat = min(mask.ylevels)
        max_lat = max(mask.ylevels)
        # Compute maximum value for degrees
        max_deg = min([(max_lon - min_lon), (max_lat - min_lat)])
        if degrees <= 0:
            raise ValueError("degrees must be greater than 0.")
        if degrees > max_deg:
            raise ValueError("The degrees value of %g is too big for this mask (maximum: %g)" % (degrees, max_deg))
        # Check if starting point is inside the original mask
        if (start_lon > max_lon) or (start_lon < min_lon):
            raise ValueError("Invalid longitude %g (min: %g, max: %g)" % (start_lon, min_lon, max_lon))
        if (start_lat > max_lat) or (start_lat < min_lat):
            raise ValueError("Invalid latitude %g (min: %g, max: %g)" % (start_lat, min_lat, max_lat))
        output = list()
        # Bottom Left point
        BL_point = [start_lon, start_lat]
        # Top Right point
        TR_point = [start_lon + degrees, start_lat + degrees]
        # Section indices
        lon_in = 0
        lat_in = 0
        # Latitude cycle
        while (TR_point[1] <= max_lat):
            # Longitude cycle
            while (TR_point[0] <= max_lon):
                # Create the Rectagle
                rect = Rectangle(BL_point[0], TR_point[0], BL_point[1], TR_point[1])
                # Create the basin
                basin = SimpleBasin("section%d%d" % (lon_in, lat_in), rect)
                # Create the SubMask and append it to output
                output.append(SubMask(basin, maskobject=mask))
                # Increment longitude
                BL_point[0] += degrees
                TR_point[0] += degrees
            # Increment latitude
            BL_point[1] += degrees
            TR_point[1] += degrees
        return output
