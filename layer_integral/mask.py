# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import netCDF4

class Mask(object):
    """
    Defines a mask from a NetCDF file
    """
    def __init__(self, filename, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon"):
        filename = str(filename)
        try:
            dset = netCDF4.Dataset(filename)
            if maskvarname in dset.variables:
                m = dset.variables[maskvarname]
                if len(m.shape) == 4:
                    self.__mask = np.array(m[0,:,:,:], dtype=np.bool)
                elif len(m.shape) == 3:
                    self.__mask = np.array(m[:,:,:], dtype=np.bool)
                else:
                    raise ValueError("Wrong shape: %s" % (m.shape,))
                self.__shape = self.__mask.shape
            else:
                raise ValueError("maskvarname '%s' not found" % (str(maskvarname),))
            if zlevelsvar in dset.variables:
                z = dset.variables[zlevelsvar]
                if len(z.shape) != 1:
                    raise ValueError("zlevelsvar must have only one dimension")
                if not z.shape[0] in self.__shape:
                    raise ValueError("cannot match %s lenght with any of %s dimensions" % (zlevelsvar, maskvarname))
                self.__zlevels = np.array(dset.variables[zlevelsvar])
            else:
                raise ValueError("zlevelsvar '%s' not found" % (str(zlevelsvar),))
            if ylevelsmatvar in dset.variables:
                self.__ylevels = np.array(dset.variables[ylevelsmatvar][:,0])
            else:
                raise ValueError("ylevelsmatvar '%s' not found" % (str(ylevelsmatvar),))
            if xlevelsmatvar in dset.variables:
                self.__xlevels = np.array(dset.variables[xlevelsmatvar][0,:])
            else:
                raise ValueError("xlevelsmatvar '%s' not found" % (str(xlevelsmatvar),))
        except:
            raise

    @property
    def mask(self):
        return self.__mask

    @property
    def zlevels(self):
        return self.__zlevels

    @property
    def shape(self):
        return self.__shape

    def convert_lon_lat_to_indices(self, lon, lat):
        """Converts longitude and latitude to the nearest indices on the mask.

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
        """
        #Input validation
        lon = float(lon)
        lat = float(lat)
        if lon > max(self.__xlevels) or lon < min(self.__xlevels):
            raise ValueError("Invalid longitude value: %f (must be between %f and %f)" % (lon, min(self.__xlevels), max(self.__xlevels)))
        if lat > max(self.__ylevels) or lat < min(self.__ylevels):
            raise ValueError("Invalid latitude value: %f (must be between %f and %f)" % (lat, min(self.__ylevels), max(self.__ylevels)))
        d_lon = np.array(self.__xlevels - lon)
        d_lon *= d_lon
        d_lat = np.array(self.__ylevels - lat)
        d_lat *= d_lat
        return np.where(d_lon == min(d_lon))[0][0] , np.where(d_lat == min(d_lat))[0][0]
