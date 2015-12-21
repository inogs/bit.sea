# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import netCDF4

class Mask(object):
    """
    Defines a mask from a NetCDF file
    """
    def __init__(self, filename, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t"):
        filename = str(filename)
        try:
            dset = netCDF4.Dataset(filename)
            if maskvarname in dset.variables:
                m = dset.variables[maskvarname]
                if len(m.shape) == 4:
                    self._mask = np.array(m[0,:,:,:], dtype=np.bool)
                elif len(m.shape) == 3:
                    self._mask = np.array(m[:,:,:], dtype=np.bool)
                else:
                    raise ValueError("Wrong shape: %s" % (m.shape,))
                self._shape = self._mask.shape
            else:
                raise ValueError("maskvarname '%s' not found" % (str(maskvarname),))
            if zlevelsvar in dset.variables:
                z = dset.variables[zlevelsvar]
                if len(z.shape) != 1:
                    raise ValueError("zlevelsvar must have only one dimension")
                if not z.shape[0] in self._shape:
                    raise ValueError("cannot match %s lenght with any of %s dimensions" % (zlevelsvar, maskvarname))
                self._zlevels = np.array(dset.variables[zlevelsvar])
            else:
                raise ValueError("zlevelsvar '%s' not found" % (str(zlevelsvar),))
            if dzvarname in dset.variables:
                self._dz = np.array(dset.variables[dzvarname][0,:,0,0])
            else:
                raise ValueError("dzvarname '%s' not found" % (str(dzvarname),))
            if ylevelsmatvar in dset.variables:
                self._ylevels = np.array(dset.variables[ylevelsmatvar][:,0])
            else:
                raise ValueError("ylevelsmatvar '%s' not found" % (str(ylevelsmatvar),))
            if xlevelsmatvar in dset.variables:
                self._xlevels = np.array(dset.variables[xlevelsmatvar][0,:])
            else:
                raise ValueError("xlevelsmatvar '%s' not found" % (str(xlevelsmatvar),))
        except:
            raise

    @property
    def mask(self):
        return self._mask

    @property
    def xlevels(self):
        return self._xlevels

    @property
    def ylevels(self):
        return self._ylevels

    @property
    def zlevels(self):
        return self._zlevels

    @property
    def dz(self):
        return self._dz

    @property
    def shape(self):
        return self._shape

    def convert_lon_lat_to_indices(self, lon, lat):
        """Converts longitude and latitude to the nearest indices on the mask.

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
        """
        #Input validation
        lon = float(lon)
        lat = float(lat)
        if lon > max(self._xlevels) or lon < min(self._xlevels):
            raise ValueError("Invalid longitude value: %f (must be between %f and %f)" % (lon, min(self._xlevels), max(self._xlevels)))
        if lat > max(self._ylevels) or lat < min(self._ylevels):
            raise ValueError("Invalid latitude value: %f (must be between %f and %f)" % (lat, min(self._ylevels), max(self._ylevels)))
        d_lon = np.array(self._xlevels - lon)
        d_lon *= d_lon
        d_lat = np.array(self._ylevels - lat)
        d_lat *= d_lat
        return np.where(d_lon == min(d_lon))[0][0] , np.where(d_lat == min(d_lat))[0][0]

    def convert_x_y_to_lon_lat(self, x, y):
        """Converts x and y coordinates to longitude and latitude.
        """
        return (self._xlevels[x], self._ylevels[y])
