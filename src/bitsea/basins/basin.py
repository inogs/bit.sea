# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
import importlib
from inspect import currentframe

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import PathPatch

from bitsea.basins.region import EmptyRegion
from bitsea.basins.region import Region


# This is the default module for the basins, i.e. the basins of this
# module do not have a prefix inside their uuids
DEFAULT_BASIN_MODULE = "bitsea.basins.V2"


class Basin(object):
    _INSTANTIATED_BASINS = {}

    def __init__(self, name, extended_name=None):
        if "." in name:
            raise ValueError(
                'Basin name cannot contain dots; received "{}"'.format(name)
            )
        if "," in name:
            raise ValueError(
                'Basin name cannot contain commas; received "{}"'.format(name)
            )

        self.name = name
        self.extended_name = extended_name
        self.region = EmptyRegion()

        # Here we save the name of the module where the basin has been
        # instantiated
        instantiation_module = currentframe().f_back.f_back
        if instantiation_module is None:
            self.__module_name = "__main__"
        else:
            self.__module_name = instantiation_module.f_globals["__name__"]

        # Now we save the current basin inside the INSTANTIATED_BASINS dict;
        # we save the basins of each module in a different sub-dictionary
        if self.__module_name not in self.__class__._INSTANTIATED_BASINS:
            self.__class__._INSTANTIATED_BASINS[self.__module_name] = {}

        module_dict = self.__class__._INSTANTIATED_BASINS[self.__module_name]
        module_dict[self.name] = self

    def __repr__(self):
        if self.extended_name is None:
            return self.name + " basin"
        else:
            return self.extended_name

    def __iter__(self):
        raise NotImplementedError

    def is_inside(self, lon, lat):
        raise NotImplementedError

    def get_uuid(self) -> str:
        """
        Return a string that uniquely identifies this basin
        """
        module_name = self.__module_name
        if module_name == DEFAULT_BASIN_MODULE:
            return self.name

        if module_name.startswith("bitsea.basins."):
            module_name = "." + module_name[len("bitsea.basins") :]
        return "{}.{}".format(module_name, self.name)

    @staticmethod
    def load_from_uuid(uuid: str):
        if uuid.startswith(".."):
            uuid = "bitsea.basins" + uuid[1:]

        if "." not in uuid:
            basin_module = DEFAULT_BASIN_MODULE
            basin_name = uuid
        else:
            basin_name = uuid.split(".")[-1]
            basin_module = uuid[: -len(basin_name) - 1]

        # We ensure that the module that contains the basin has been loaded
        if basin_module != "__main__":
            importlib.import_module(basin_module)

        return Basin._INSTANTIATED_BASINS[basin_module][basin_name]

    def cross(self, region_or_basin):
        if isinstance(region_or_basin, Basin):
            return self.region.cross(region_or_basin.region)
        elif isinstance(region_or_basin, Region):
            return self.region.cross(region_or_basin)

    def grid(
        self,
        lon_window=(-8, 38),
        lat_window=(30, 47),
        lon_points=100,
        lat_points=100,
    ):
        """
        Generate a grid of equally spaced points inside a rectangular domain
        and check which ones of those points are inside the current basin
        """
        lon_grid = np.linspace(
            lon_window[0], lon_window[1], lon_points, endpoint=True
        )

        lat_grid = np.linspace(
            lat_window[0], lat_window[1], lat_points, endpoint=True
        )

        inside_domain = self.is_inside(
            lon_grid, lat_grid.reshape((lat_points, 1))
        )
        return lon_grid, lat_grid, inside_domain

    def plot(
        self,
        lon_window=(-8, 38),
        lat_window=(30, 47),
        lon_points=100,
        lat_points=100,
        color="tab:blue",
        alpha=1,
        zorder=None,
        axis=None,
        fill=True,
        transform=None,
    ):
        lon_grid, lat_grid, inside_domain = self.grid(
            lon_window, lat_window, lon_points, lat_points
        )

        if axis is None:
            axis = plt.gca()

        kwargs = {}

        if fill:
            plot_f = axis.contourf
            cm = LinearSegmentedColormap.from_list(
                "none", ["black", color], N=2
            )
            kwargs["cmap"] = cm
        else:
            plot_f = axis.contour
            kwargs["colors"] = [color]

        if transform is not None:
            kwargs["transform"] = transform

        if zorder is not None:
            kwargs["zorder"] = zorder

        current_plot = plot_f(
            lon_grid, lat_grid, inside_domain, [0.1, 1], alpha=alpha, **kwargs
        )

        return current_plot


class SimpleBasin(Basin):
    def __init__(self, name, region, extended_name=None):
        super(SimpleBasin, self).__init__(name, extended_name)
        self.region = region

    def __iter__(self):
        return [self].__iter__()

    def is_inside(self, lon, lat):
        return self.region.is_inside(lon=lon, lat=lat)


class SimplePolygonalBasin(SimpleBasin):
    @property
    def borders(self):
        return self.region.borders

    def plot(
        self,
        lon_window=(-8, 38),
        lat_window=(30, 47),
        lon_points=100,
        lat_points=100,
        color="tab:blue",
        alpha=1,
        zorder=None,
        axis=None,
        fill=True,
        transform=None,
    ):
        patch_kwargs = {"color": color, "alpha": alpha, "fill": fill}
        if zorder is not None:
            patch_kwargs["zorder"] = zorder
        if transform is not None:
            patch_kwargs["transform"] = transform

        patch = PathPatch(self.region.path, **patch_kwargs)

        if axis is None:
            axis = plt.gca()

        axis.add_patch(patch)

        return patch


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

    def plot(
        self,
        lon_window=(-8, 38),
        lat_window=(30, 47),
        lon_points=100,
        lat_points=100,
        color="tab:blue",
        alpha=1,
        zorder=None,
        axis=None,
        fill=True,
        transform=None,
    ):
        for basin in self.basin_list:
            basin.plot(
                lon_window=lon_window,
                lat_window=lat_window,
                lon_points=lon_points,
                lat_points=lat_points,
                color=color,
                alpha=alpha,
                zorder=zorder,
                axis=axis,
                fill=fill,
                transform=transform,
            )

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
