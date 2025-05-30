import csv
import pathlib
import re
import sys
from abc import ABC
from abc import abstractmethod
from os import PathLike
from typing import Dict
from typing import TextIO
from typing import Union

import numpy as np
from matplotlib.path import Path

from bitsea.utilities.string_manipulation import ComposedString
from bitsea.utilities.string_manipulation import read_csv_txt


class Region(ABC):
    @abstractmethod
    def is_inside(self, *, lon, lat):
        raise NotImplementedError

    def __add__(self, r):
        return RegionUnion(self, r)

    def intersect(self, r):
        return RegionIntersection(self, r)


class EmptyRegion(Region):
    def is_inside(self, *, lon, lat):
        if hasattr(lon, "__len__") or hasattr(lat, "__len__"):
            lon, lat = np.broadcast_arrays(np.asarray(lon), np.asarray(lat))
            return np.zeros_like(lon, dtype=bool)
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

        # Compute the min and the max of each variable
        self.__lon_range = (min(lon_list), max(lon_list))
        self.__lat_range = (min(lat_list), max(lat_list))

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

        lon_in_range = np.logical_and(
            points_coord[:, 0] >= self.__lon_range[0],
            points_coord[:, 0] <= self.__lon_range[1],
        )
        lat_in_range = np.logical_and(
            points_coord[:, 1] >= self.__lat_range[0],
            points_coord[:, 1] <= self.__lat_range[1],
        )

        inside_rectangle = np.logical_and(lon_in_range, lat_in_range)

        if np.any(inside_rectangle):
            inside = np.zeros_like(inside_rectangle)
            inside[inside_rectangle] = self.path.contains_points(
                points_coord[inside_rectangle, :]
            )
        else:
            inside = inside_rectangle

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
        return tuple(p for p in zip(self.__lon_list, self.__lat_list))

    def cross(self, another_region):
        # The following lines are useful if another_region is a basin
        if hasattr(another_region, "region"):
            another_region = another_region.region

        if isinstance(another_region, Polygon):
            return self.path.intersects_path(another_region.path, filled=True)
        elif isinstance(another_region, RegionUnion):
            return another_region.cross(self)
        elif isinstance(another_region, EmptyRegion):
            return False

        return NotImplemented

    @staticmethod
    def from_WKT_file(filepath: Union[PathLike, str]):
        # Ensure filepath is a Path and not a simple string
        filepath = pathlib.Path(filepath)

        poly_list = []
        csv.field_size_limit(sys.maxsize)

        with filepath.open("r") as f:
            reader = csv.reader(f, delimiter="\t")
            for line in reader:
                poly_list.append(line[0])

        string_poly = " ".join(poly_list)
        string_poly_sel = string_poly[15:]

        single_polygons = re.findall(r"\((.*?)\)", string_poly_sel)
        n_poly = len(single_polygons)
        if n_poly > 1:
            # case to be implemented
            raise ValueError("filename containing more than one polygon")

        pat = re.compile(r"(-*\d+\.\d+ -*\d+\.\d+),*")
        matches = pat.findall(single_polygons[0])

        if matches is None:
            raise ValueError("Invalid string inside the polygon file")

        lst = [tuple(map(float, m.split())) for m in matches]
        n_coords = len(lst)
        lon = []
        lat = []

        for ii in range(n_coords):
            lon.append(lst[ii][0])
            lat.append(lst[ii][1])
        return Polygon(lon, lat)

    @staticmethod
    def read_WKT_polygon(
        polygon_description: str, allowed_classes=("POLYGON", "LINESTRING")
    ):
        """
        Converts a polygon description string (which starts with "POLYGON" and
        contains the coordinates of the polygon) into an instance of this
        class.
        """
        if "(" not in polygon_description:
            raise ValueError(
                'The polygon description must start with "POLYGON" (or with '
                "the name of the shape class) and contain the coordinates of "
                "the polygon starting with a ("
            )
        parenthesis_start = polygon_description.find("(")
        polygon_class = polygon_description[:parenthesis_start].strip()
        if polygon_class.lower() not in (c.lower() for c in allowed_classes):
            raise ValueError(
                "The polygon description must start with the name of the "
                f"polygon class, which may be one of {allowed_classes}; "
                f"received {polygon_class}"
            )
        point_str = ComposedString.split_on_parenthesis(
            polygon_description[parenthesis_start:].strip()
        )
        # Possible errors to handle:
        # Additional text exists after the first POLYGON(...) definition.
        if len(point_str.elements) > 1:
            raise ValueError(
                "Invalid string after the first parenthesis: "
                f"{polygon_description}"
            )
        # The content inside POLYGON( ) does not contain additional parentheses.
        if not isinstance(point_str.elements[0], ComposedString):
            raise ValueError(
                "The first argument of the polygon description must be a "
                f"tuple of points; received {polygon_description}"
            )
        # The content inside POLYGON() contains multiple nested parentheses.
        if len(point_str.elements[0].elements) != 1:
            raise ValueError(
                "The arguments list of the polygon description must be "
                f"one single tuple of points; received {polygon_description}"
            )

        # Extract a list of points from the parentheses.
        # Each point is in the format: "lon1 lat1, lon2 lat2, lon3 lat3, ..."
        point_list = str(point_str.elements[0].elements[0])
        points = point_list.split(",")
        lon = []
        lat = []
        for point in points:
            point_split = point.split()
            lon.append(float(point_split[0]))
            lat.append(float(point_split[1]))
        return Polygon(lon, lat)

    @staticmethod
    def read_WKT_file(file_pointer: TextIO) -> Dict:
        """
        This method reads a WKT file (like the ones that can be downloaded from
        Google My Maps) and returns a dictionary that associates to each
        polygon name its corresponding polygon.

        This function is really similar to the `from_WKT_file method`, but it
        can read more than one polygon. On the other hand, the polygon that
        this function reads are simple polygons and not Multipolygons.
        """
        csv_data = read_csv_txt(file_pointer.read())

        if len(csv_data) == 0:
            return {}

        if len(csv_data[0]) != 3:
            raise ValueError(
                "The WKT file must contain three columns: "
                "polygon coordinates, polygon name and polygon description"
            )

        polys = {}
        for i, line in enumerate(csv_data):
            # Skip the first line if contains the header
            if i == 0 and line[0].lower() == "wkt":
                continue
            poly_name = line[1]
            poly_raw = line[0]
            poly = Polygon.read_WKT_polygon(poly_raw)
            polys[poly_name] = poly

        return polys


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
    This region is defined as a 2D polygon on the surface coupled with a
    bathymetric condition. The region is the part of the polygon where the
    bathymetry satisfies the condition.

    By definition, lower_than means "all the points whose bathymetry is STRICTLY
    greater than submitted number", while upper_than means "all the points whose
    bathymetry is smaller OR EQUAL than the submitted number". This avoids
    intersections between basins when using the same pivot numbers (for example,
    a basin with lower_than=100 and another with upper_than=100).
    """

    def __init__(
        self,
        lon_list,
        lat_list,
        bathymetry,
        deeper_than=None,
        shallower_than=None,
    ):
        self.polygon = Polygon(lon_list, lat_list)

        self.bathymetry = bathymetry

        if deeper_than is None and shallower_than is None:
            raise ValueError(
                "At least one between bathymetric_min and bathymetric_max"
                "must be different from None. Otherwise, simply use a Polygon"
            )

        self.deeper_than = deeper_than
        self.shallower_than = shallower_than

    def is_inside(self, lon, lat):
        inside_poly = self.polygon.is_inside(lon, lat)

        point_depth = self.bathymetry(lon, lat)

        bathymetry_domain = self.bathymetry.is_inside_domain(lon, lat)

        if self.deeper_than is not None:
            bathymetric_ok_min = point_depth > self.deeper_than
        else:
            bathymetric_ok_min = True

        if self.shallower_than is not None:
            bathymetric_ok_max = point_depth <= self.shallower_than
        else:
            bathymetric_ok_max = True

        bathymetric_ok = np.logical_and(bathymetric_ok_min, bathymetric_ok_max)
        bathymetric_ok = np.logical_and(bathymetric_ok, bathymetry_domain)

        return np.logical_and(inside_poly, bathymetric_ok)


class RegionUnion(Region):
    def __init__(self, r1, r2):
        self._r1 = r1
        self._r2 = r2

    def is_inside(self, lon, lat):
        ins1 = self._r1.is_inside(lon, lat)
        ins2 = self._r2.is_inside(lon, lat)
        return np.logical_or(ins1, ins2)

    def cross(self, another_region):
        cross1 = self._r1.cross(another_region)
        cross2 = self._r2.cross(another_region)
        return bool(cross1 + cross2)


class RegionIntersection(Region):
    def __init__(self, *args):
        self._regions = tuple(args)

    def is_inside(self, lon, lat):
        if len(self._regions) == 0:
            return np.full_like(np.broadcast(lon, lat), True, dtype=bool)

        lon, lat = np.broadcast_arrays(lon, lat)

        # This is the case when lon and lat are just numbers and not an array
        if len(lon.shape) == 0:
            out_value = self._regions[0].is_inside(lon, lat)
            for r in self._regions[1:]:
                if not out_value:
                    break
                out_value = r.is_inside(lon, lat)
            return out_value

        out_values = self._regions[0].is_inside(lon, lat)
        for r in self._regions[1:]:
            out_values[out_values] = r.is_inside(
                lon[out_values], lat[out_values]
            )
        return out_values
