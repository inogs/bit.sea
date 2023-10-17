from collections.abc import Callable
from abc import ABC, abstractmethod

from netCDF4 import Dataset
import numpy as np
from scipy import interpolate


class Bathymetry(ABC, Callable):
    """
    A Bathymetry is a callable object that returns a bathymetry for a specific
    point (described by its longitude and its latitude). A call to a Bathymetry
    object always returns a non-negative number.

    Each Bathymetry has an associated domain that represents the geographical
    region where the object knows the bathymetry. When a point is outside the
    domain, the bathymetry returns always 0. If you want to check if a point is
    inside the domain, you can use the method `is_inside_domain`.

    Both the call method of the bathymetry and the `is_inside_domain` method
    must support the numpy vectorization and the broadcasting.
    """
    @abstractmethod
    def __call__(self, lon, lat):
        raise NotImplementedError

    def is_inside_domain(self, lon, lat):
        raise NotImplementedError


class GEBCOBathymetry(Bathymetry):
    def __init__(self, data_file_path):
        self._path = data_file_path
        with Dataset(self._path, 'r') as f:
            self._lon = np.array(f.variables['lon'][:], dtype=np.float64)
            self._lat = np.array(f.variables['lat'][:], dtype=np.float64)
            self._values = np.array(
                f.variables['elevation'][:],
                dtype=np.float32
            )

            # Get the masked values
            values_mask = np.ma.getmask(f.variables['elevation'][:])

        # Set the masked values to 0
        self._values[values_mask] = 0

        # Now we set the bathymetry to 0 on the ground and we reverse the sign
        self._values[self._values > 0] = 0
        self._values *= -1

        self._interpolator = interpolate.RegularGridInterpolator(
            (self._lat, self._lon),
            self._values,
            bounds_error=False,
            fill_value=0
        )

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, repr(self._path))

    def __call__(self, lon, lat):
        return self._interpolator((lat, lon))

    def is_inside_domain(self, lon, lat):
        inside_lon = np.logical_and(
            lon >= self._lon[0],
            lon <= self._lon[-1],
        )
        inside_lat = np.logical_and(
            lat >= self._lat[0],
            lat <= self._lat[-1],
        )
        return np.logical_and(inside_lon, inside_lat)
