# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
import numpy as np
from matplotlib.path import Path


class Region(object):
    def is_inside(self, lon, lat):
        raise NotImplementedError
    
    def __add__(self, r):
        return RegionUnion(self, r)
    

class EmptyRegion(Region):
    def is_inside(self, lon, lat):
        if hasattr(lon, "__len__"):
            assert len(lon) == len(lat)
            return np.zeros((len(lon),), dtype=np.bool_)
        else:
            return False

    def cross(self, another_region):
        return False


class Polygon(Region):
    def __init__(self, lon_list, lat_list):
        assert len(lon_list) == len(lat_list)
        assert len(lon_list) > 2

        # Check that input is a list
        lon_list = list(lon_list)
        lat_list = list(lat_list)
        
        # Save the longitude and the latitude lists
        self.__lon_list = lon_list
        self.__lat_list = lat_list

        # Ensure that the input is close
        if lon_list[-1] != lon_list[0] or lat_list[-1] != lat_list[0]:
            lon_list.append(lon_list[0])
            lat_list.append(lat_list[0])
        
        # Create a path object
        codes = [Path.LINETO] * len(lon_list)
        codes[0] = Path.MOVETO
        codes[-1] = Path.CLOSEPOLY
        
        coords = [(lon_list[i], lat_list[i]) for i in range(len(lon_list))]
        
        self.path = Path(coords, codes)
        
    def is_inside(self, lon, lat):
        points_coord = np.stack(np.broadcast_arrays(lon, lat), axis=-1)
        assert points_coord.shape[-1] == 2

        # path.contains_points accepts only 1d arrays of pairs. reshaping for
        # single numbers and more complicated shapes
        reshaped = False
        undo_reshaping = None
        # Check if we have a vector of just two elements (i.e., shape (2,))
        if len(points_coord.shape) < 2:
            points_coord = points_coord.reshape(1, 2)
            reshaped = True

            def undo_reshaping(x):
                return x[0]

        elif len(points_coord.shape) > 2:
            current_shape = points_coord.shape
            points_coord = points_coord.reshape((-1, 2))
            reshaped = True

            def undo_reshaping(x):
                return x.reshape(current_shape[:-1])

        inside = self.path.contains_points(points_coord)

        if reshaped:
            return undo_reshaping(inside)
        else:
            return inside
    
    @property 
    def border_latitudes(self):
        return self.__lat_list

    @property 
    def border_longitudes(self):
        return self.__lon_list

    @property 
    def borders(self):
        output = []
        for i in range(len(self.__lon_list)):
            lon = self.border_longitudes[i]
            lat = self.border_latitudes[i]
            output.append((lon, lat))
        return tuple(output)

    def cross(self, another_region):
        # The following lines are useful if another_region is a basin
        if hasattr(another_region, "region"):
            another_region = another_region.region

        if isinstance(another_region, Polygon):
            return np.bool_(self.path.intersects_path(another_region.path, filled=True))
        elif isinstance(another_region, RegionUnion):
            return another_region.cross(self)
        elif isinstance(another_region, EmptyRegion):
            return False

        return NotImplemented


class Rectangle(Polygon):
    def __init__(self, lonmin, lonmax, latmin, latmax):
        self.lonmin = lonmin
        self.lonmax = lonmax
        self.latmin = latmin
        self.latmax = latmax

        lonlist = [lonmin, lonmax, lonmax, lonmin, lonmin]
        latlist = [latmin, latmin, latmax, latmax, latmin]

        super(Rectangle, self).__init__(lonlist, latlist)

    def is_inside(self, lon, lat):
        lat_inside = np.logical_and(lat >= self.latmin, lat <= self.latmax)
        lon_inside = np.logical_and(lon >= self.lonmin, lon <= self.lonmax)
        inside = np.logical_and(lat_inside, lon_inside)
        return inside


class BathymetricPolygon(Region):
    """
    This region is defined as a 2D polygon on the surface and with a bathymetric
    condition. The region is the part of the polygon where the bathymetry
    satisfies the condition
    """
    def __init__(self, lon_list, lat_list, bathymetric_min=None,
                 bathymetric_max=None):
        self.polygon = Polygon(lon_list, lat_list)

        if bathymetric_min is None and bathymetric_max is None:
            raise ValueError(
                'At least one between bathymetric_min and bathymetric_max must '
                'be different from None. Otherwise, simply use a Polygon'
            )

        self.bathymetric_min = bathymetric_min
        self.bathymetric_max = bathymetric_max

    def is_inside(self, lon, lat):
        inside_poly = self.polygon.is_inside(lon, lat)

        # TODO: find a way to compute points depth
        # Possible optimization: compute the point_depth only for the points
        # that are inside_poly and set a 0 for the others
        point_depth = np.ones_like(lon)

        if self.bathymetric_min is not None:
            bathymetric_ok_min = point_depth > self.bathymetric_min
        else:
            bathymetric_ok_min = np.ones_like(lon, dtype=bool)

        if self.bathymetric_max is not None:
            bathymetric_ok_max = point_depth < self.bathymetric_max
        else:
            bathymetric_ok_max = np.ones_like(lon, dtype=bool)

        bathymetric_ok = np.logical_and(
            bathymetric_ok_min,
            bathymetric_ok_max
        )

        return np.logical_and(inside_poly, bathymetric_ok)


class RegionUnion(Region):
    def __init__(self, r1, r2):
        self._r1 = r1
        self._r2 = r2
    
    def is_inside(self, lon, lat):
        ins1 = self._r1.is_inside(lon,lat)
        ins2 = self._r2.is_inside(lon,lat)
        return np.bool_(ins1 + ins2)

    def cross(self, another_region):
        cross1 = self._r1.cross(another_region)
        cross2 = self._r2.cross(another_region)
        return np.bool_(cross1+cross2)
