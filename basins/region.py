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
        if hasattr(lon,"__len__"):
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
        codes[0]  = Path.MOVETO
        codes[-1] = Path.CLOSEPOLY
        
        coords = [ (lon_list[i], lat_list[i]) for i in range(len(lon_list))]
        
        self.path = Path(coords, codes)
        
    def is_inside(self, lon, lat):
        points_coord = np.array((lon,lat)).T

        reshaped = False
        if len(points_coord.shape) < 2:
            points_coord = points_coord.reshape(1,2)
            reshaped = True

        inside = self.path.contains_points(points_coord)

        if reshaped:
            return inside[0]
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
        else:
            raise NotImplementedError


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
        lat_inside = (lat>= self.latmin) * (lat <= self.latmax) 
        lon_inside = (lon >= self.lonmin) * (lon <= self.lonmax)
        inside = np.bool_(lat_inside * lon_inside)
        return inside



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
