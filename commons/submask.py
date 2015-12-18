# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import netCDF4
from commons.mask import *
from basins.basin import Basin

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
            super(SubMask, self).__init__(filename, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t")
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
