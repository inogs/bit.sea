import numpy as np
from matplotlib.path import Path

class Region():
    def is_inside(self, lon, lat):
        raise NotImplementedError
    
    def __add__(self, r):
        return RegionUnion(self, r)
    
    def __mul__(self, r):
        return RegionIntersection(self,r)
    
    def __neg__(self):
        return RegionComplement(self)
    
    def __sub__(self, r):
        return self * (-r)


class Rectangle(Region):
    def __init__(self, lonmin, lonmax, latmin, latmax):
        self.lonmin = lonmin
        self.lonmax = lonmax
        self.latmin = latmin
        self.latmax = latmax
    
    def is_inside(self, lon, lat):
        lat_inside = (lat>= self.latmin) * (lat <= self.latmax) 
        lon_inside = (lon >= self.lonmin) * (lon <= self.lonmax)
        inside = np.bool_(lat_inside * lon_inside)
        return inside


class Circle(Region):
    def __init__(self, center_lon, center_lat, radious):
        self.center_lon = center_lon
        self.center_lat = center_lat
        self.radious = radious
    
    def is_inside(self, lon, lat):
        distance = (lon - self.center_lon)**2 + (lat - self.center_lat)**2
        return distance < self.radious**2


class Polygon(Region):
    def __init__(self, lon_list, lat_list):
        assert len(lon_list) == len(lat_list)
        assert len(lon_list) > 2

        # Check that input is a list
        lon_list = list(lon_list)
        lat_list = list(lat_list)
        
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


class RegionUnion(Region):
    def __init__(self, r1, r2):
        self._r1 = r1
        self._r2 = r2
    
    def is_inside(self, lon, lat):
        ins1 = self._r1.is_inside(lon,lat)
        ins2 = self._r2.is_inside(lon,lat)
        return np.bool_(ins1 + ins2)


class RegionIntersection(Region):
    def __init__(self, r1, r2):
        self._r1 = r1
        self._r2 = r2
    
    def is_inside(self, lon, lat):
        ins1 = self._r1.is_inside(lon,lat)
        ins2 = self._r2.is_inside(lon,lat)
        return np.bool_(ins1 * ins2)


class RegionComplement(Region):
    def __init__(self, r):
        self._r = r

    def is_inside(self, lon, lat):
        try:
            inside = not self._r.is_inside(lon, lat)
        except ValueError:
            inside = ~self._r.is_inside(lon, lat)
        
        return inside
