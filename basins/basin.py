# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

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

    def grid(self, lon_window=(-8, 38), lat_window=(30, 47), lon_points=100,
             lat_points=100):
        """
        Generate a grid of equally spaced points inside a rectangular domain
        and check which ones of those points are inside the current basin
        """
        lon_grid = np.linspace(
            lon_window[0],
            lon_window[1],
            lon_points,
            endpoint=True
        )

        lat_grid = np.linspace(
            lat_window[0],
            lat_window[1],
            lat_points,
            endpoint=True
        )

        inside_domain = self.is_inside(
            lon_grid,
            lat_grid.reshape((lat_points, 1))
        )
        return lon_grid, lat_grid, inside_domain

    def plot(self, lon_window=(-8, 38), lat_window=(30, 47), lon_points=100,
             lat_points=100, color="tab:blue", alpha=1, axes=None, filled=True):
        lon_grid, lat_grid, inside_domain = self.grid(
            lon_window,
            lat_window,
            lon_points,
            lat_points
        )

        if axes is None:
            axes = plt.gca()

        kwargs = {}

        if filled:
            plot_f = axes.contourf
            cm = LinearSegmentedColormap.from_list(
                "none",
                ['black', color],
                N=2
            )
            kwargs['cmap'] = cm
        else:
            plot_f = axes.contour
            kwargs['colors'] = [color]

        current_plot = plot_f(
            lon_grid,
            lat_grid,
            inside_domain,
            [0.1, 1],
            alpha=alpha,
            **kwargs
        )

        return current_plot



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


class SimpleBathymetricBasin(SimpleBasin):
    pass


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
            if hasattr(lon, "__len__") or hasattr(lat, "__len__"):
                output_shape = np.broadcast(lon, lat).shape
                return np.zeros(shape=output_shape, dtype=bool)
            else:
                return False

        output = self.basin_list[0].is_inside(lon, lat)
        for current_basin in self.basin_list[1:]:
            output = np.logical_or(output, current_basin.is_inside(lon, lat))
        return output
