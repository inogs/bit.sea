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
                self._ylevels = np.array(dset.variables[ylevelsmatvar])
            else:
                raise ValueError("ylevelsmatvar '%s' not found" % (str(ylevelsmatvar),))
            if xlevelsmatvar in dset.variables:
                self._xlevels = np.array(dset.variables[xlevelsmatvar])
            else:
                raise ValueError("xlevelsmatvar '%s' not found" % (str(xlevelsmatvar),))
            e1t = np.array(dset.variables['e1t'][0,0,:,:]).astype(np.float32)
            e2t = np.array(dset.variables['e2t'][0,0,:,:]).astype(np.float32)
            self._area = e1t*e2t
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
    @property
    def area(self):
        return self._area

    def convert_lon_lat_to_indices(self, lon, lat):
        """Converts longitude and latitude to the nearest indices on the mask.

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
        Returns: a tuple of numbers, the first one is the longitude index and
        the other one is the latitude index.
        """
        #Input validation
        lon = float(lon)
        lat = float(lat)
        min_lon = self._xlevels.min()
        max_lon = self._xlevels.max()
        min_lat = self._ylevels.min()
        max_lat = self._ylevels.max()
        if lon > max_lon or lon < min_lon:
            raise ValueError("Invalid longitude value: %f (must be between %f and %f)" % (lon, min_lon, max_lon))
        if lat > max_lat or lat < min_lat:
            raise ValueError("Invalid latitude value: %f (must be between %f and %f)" % (lat, min_lat, max_lat))
        #Longitude distances matrix
        d_lon = np.array(self._xlevels - lon)
        d_lon *= d_lon
        #Latitude distances matrix
        d_lat = np.array(self._ylevels - lat)
        d_lat *= d_lat
        #Compute minimum indices
        min_d_lon = d_lon.min()
        min_d_lat = d_lat.min()
        lon_indices = np.where(d_lon == min_d_lon)
        lat_indices = np.where(d_lat == min_d_lat)
        return lon_indices[1][0], lat_indices[0][0]

    def convert_x_y_to_lon_lat(self, x, y):
        """Converts x and y coordinates to longitude and latitude.
        """
        return (self._xlevels[y,x], self._ylevels[y,x])

    def getDepthIndex(self, z):
        '''Converts a depth expressed in meters in the corresponding index level
        The returned value is an integer indicating the previous (not nearest) depth in the z-levels.

        Example:
        M = Mask(filename)
        k = M.getDepthIndex(200.)
        M.zlevels[k]
          returns 192.60
        '''
        jk_m = 0
        for jk,depth in enumerate(self.zlevels):
            if depth < z:
                jk_m=jk
        return jk_m

    def mask_at_level(self,z):
        '''
        Returns a 2d map of logicals, for the depth (m) provided as argument.
        as a slice of the tmask.
        '''
        jk_m = self.getDepthIndex(z)
        level_mask = self.mask[jk_m,:,:].copy()
        return level_mask
