# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
import numpy as np
from basins.region import Region, EmptyRegion

class Basin(object):

    def __init__(self, name, extended_name=None):
        self.name = name
        self.extended_name = extended_name
        self.region = EmptyRegion()

    def __repr__(self):
        if self.extended_name is None:
            return self.name + ' basin'
        else:
            return self.extended_name

    def __iter__(self):
        raise NotImplementedError

    def is_inside(self, lon, lat):
        raise NotImplementedError

    def cross(self, region_or_basin):
        if isinstance(region_or_basin, Basin):
            return self.region.cross(region_or_basin.region)
        elif isinstance(region_or_basin, Region):
            return self.region.cross(region_or_basin)


class SimpleBasin(Basin):
    def __init__(self, name, region, extended_name=None):
        super(SimpleBasin, self).__init__(name, extended_name)
        self.region = region
    
    def __iter__(self):
        return [self].__iter__()

    def is_inside(self, lon, lat):
        return self.region.is_inside(lon, lat)

class SimplePolygonalBasin(SimpleBasin):
    @property
    def borders(self):
        return self.region.borders


class ComposedBasin(Basin):
    def __init__(self, name, basin_list, extended_name=None):
        super(ComposedBasin, self).__init__(name, extended_name)
        self.basin_list = basin_list
        if len(basin_list) > 0:
            region = basin_list[0].region
            for i in range(1, len(basin_list)):
                region = region + basin_list[i].region
            self.region = region
    
    def __iter__(self):
        return self.basin_list.__iter__()
    
    def is_inside(self, lon, lat):
        if len(self.basin_list) == 0:
            if hasattr(lon,"__len__"):
                assert len(lon) == len(lat)
                return np.zeros((len(lon),), dtype=np.bool_)
            else:
                return False
        else:
            output = self.basin_list[0].is_inside(lon, lat)
            for i in range(1, len(self.basin_list)):
                output = np.logical_or(output, self.basin_list[i].is_inside(lon, lat))
            return output
