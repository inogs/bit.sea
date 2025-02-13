from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from enum import Enum
from itertools import product as cart_prod
from typing import Optional
from typing import Tuple

import numpy as np
from geographiclib.geodesic import Geodesic
from numpy.typing import ArrayLike


class GridSidesAlgorithm(Enum):
    NEMO = 1
    GREAT_CIRCLE = 2
    GEODESIC = 3


@dataclass
class Ellipsoid:
    """
    An Ellipsoid defines the shape of the Earth; if `a` is the major axis and
    `b` is the minor axis of the ellipsoid, the flattening is defined as
    `(a - b) / a`
    """

    name: str
    major_axis: float
    flattening: float

    @property
    def minor_axis(self):
        return (1 - self.flattening) * self.major_axis


# This is the standard WGS84 ellipsoid (cfr.
# https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84)
WGS84 = Ellipsoid("WGS84", 6378137, 1 / 298.257223563)

# Define an object from the geographiclib library that can solve
# the geodesic problem on our ellipsoid
GEODESIC = Geodesic(WGS84.major_axis, WGS84.flattening)


def _geodetic_distance(
    p1: Tuple[float, float], p2: Tuple[float, float]
) -> float:
    """
    Calculate the geodetic distance between two points on Earth.

    This function computes the geodetic distance between two geographical
    points, specified by their latitude and longitude, using the Geodesic
    Inverse method.
    It accurately calculates the shortest distance over the Earth's
    ellipsoidal surface.

    This function returns the distance in meters.

    This function is not intended to be called directly. Instead, it is
    used by the `GeoPySidesCalculator` class.

    Args:
        p1: The latitude and longitude of the first point in degrees.
        p2: The latitude and longitude of the second point in degrees.

    Returns:
        float: The geodetic distance between the two specified points,
            in meters.
    """
    lat1, lon1 = p1
    lat2, lon2 = p2
    return GEODESIC.Inverse(lat1, lon1, lat2, lon2, Geodesic.DISTANCE)["s12"]


def _great_circle_distance(
    p1: Tuple[float, float], p2: Tuple[float, float]
) -> float:
    """
    Calculate the great-circle distance between two points on the Earth's
    surface.

    This function computes the great-circle distance, which is the shortest
    distance between two points on the surface of a sphere, in this case,
    the Earth. The calculation is performed using the haversine formula.
    The Earth is approximated as a sphere of radius 6,371,009 meters.

    This function returns the distance in meters.

    This function is not intended to be called directly. Instead, it is
    used by the `GeoPySidesCalculator` class.

    Args:
        p1: The latitude and longitude of the first point in degrees.
        p2: The latitude and longitude of the second point in degrees.

    Returns:
        float: The great-circle distance between the two points, in meters.
    """
    lat1, lon1 = np.radians(p1[0]), np.radians(p1[1])
    lat2, lon2 = np.radians(p2[0]), np.radians(p2[1])

    sin_lat1, cos_lat1 = np.sin(lat1), np.cos(lat1)
    sin_lat2, cos_lat2 = np.sin(lat2), np.cos(lat2)

    delta_lon = lon2 - lon1
    cos_delta_lon, sin_delta_lon = np.cos(delta_lon), np.sin(delta_lon)

    d = np.arctan2(
        np.sqrt(
            (cos_lat2 * sin_delta_lon) ** 2
            + (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_delta_lon) ** 2
        ),
        sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lon,
    )

    return 6371009 * d


def extend_from_average(
    v: np.ndarray, axis: int = 0, first_value: Optional[ArrayLike] = None
) -> np.ndarray:
    """
    Given a vector `v`, create a new vector `w` with the same shape as
    `v` on all axes except for the specified `axis`, where `w` has one
    additional element. The vector `w`
    is constructed such that:

        v[i] = (w[i] + w[i + 1]) / 2

    for each index `i` along the specified axis.

    This is especially useful for calculating the positions of the
    faces of a grid based on the positions of its centers.

    Args:
        v (np.ndarray): the array with the positions
        axis (int): the axis along which to extend the vector
        first_value (Optional[ArrayLike]): the first value of the output vector;
          if this is `None`, then the first value is chosen such that w[1] is
          in the middle between v[0] and v[1].
    """
    if axis < 0:
        axis += v.ndim

    if axis >= v.ndim:
        raise IndexError(
            f"Invalid axis {axis} when the array has only {v.ndim} dimensions."
        )
    if axis < 0:
        raise IndexError(
            f"Invalid axis {axis - v.ndim} when the array has only {v.ndim} "
            f"dimensions."
        )

    v = v.view()
    if axis != 0:
        v = np.moveaxis(v, axis, 0)

    output = np.empty((v.shape[0] + 1,) + v.shape[1:], dtype=v.dtype)

    if first_value is None:
        output[1] = (v[0] + v[1]) / 2.0
        output[0] = 2 * v[0] - output[1]
    else:
        output[0] = first_value

    for i in range(1, v.shape[0] + 1):
        output[i] = 2 * v[i - 1] - output[i - 1]

    if axis != 0:
        output = np.moveaxis(output, 0, axis)

    return output


class SidesCalculator(ABC):
    """A `SideCalculator` computes the size of the grid cell sides based on the
    coordinates of their centers. It provides the size of the longitudinal side
    (east-west) and the latitudinal side (north-south) for each cell in the
    grid.

    Since having only the positions of the cell centers is insufficient
    to precisely reconstruct the grid's shape, the concept of "size of
    the longitudinal side of a cell" is an approximation (as the length
    can vary across the cell). To accommodate this, various
    implementations are provided, each making different assumptions and
    employing different algorithms.

    All the distances returned by a SidesCalculator are in meters.
    """

    @abstractmethod
    def __call__(
        self, xlevels: np.ndarray, ylevels: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute the size of the cells in a grid.

        Args:
            xlevels (np.ndarray): A 2D array containing the longitudes
              of the grid cell centers.
            ylevels (np.ndarray): A 2D array containing the latitudes
              of the grid cell centers.

        Returns:
            tuple: Two numpy arrays, `e1t` and `e2t`, where:
              - `e1t[i, j]` is the length of the segment of the
                parallel passing through the cell centered at
                `(xlevels[i, j], ylevels[i, j])` that lies within the
                cell.
              - `e2t[i, j]` is the corresponding length along the
                meridian.
        """
        raise NotImplementedError


class NemoGridSidesCalculator(SidesCalculator):
    """This object is a reimplementation of the algorithm used by the NEMO
    model to compute the size of grid cells along longitude and latitude when
    the grid is regular.

    It corresponds to "case 1" in the following file from the NEMO
    source code:
    `domhgr.F90 <https://forge.ipsl.jussieu.fr/nemo/browser/utils/tools/DOMAINcfg/src/domhgr.F90?rev=13204>`_.
    """

    EARTH_RADIUS = 6371229

    def __call__(self, xlevels, ylevels) -> Tuple[np.ndarray, np.ndarray]:
        cell_lon_size = np.empty_like(xlevels, dtype=np.float32)
        cell_lat_size = np.empty_like(xlevels, dtype=np.float32)

        cell_lon_size[:, 1:-1] = (xlevels[:, 2:] - xlevels[:, :-2]) / 2
        cell_lon_size[:, 0] = xlevels[:, 1] - xlevels[:, 0]
        cell_lon_size[:, -1] = xlevels[:, -1] - xlevels[:, -2]

        cell_lat_size[1:-1, :] = (ylevels[2:, :] - ylevels[:-2, :]) / 2
        cell_lat_size[0, :] = ylevels[1, :] - ylevels[0, :]
        cell_lat_size[-1, :] = ylevels[-1, :] - ylevels[-2, :]

        lon_size_rads = np.radians(cell_lon_size)
        lat_size_rads = np.radians(cell_lat_size)
        r = self.EARTH_RADIUS
        e1t = r * np.cos(np.radians(ylevels)) * lon_size_rads
        e2t = r * lat_size_rads

        return e1t, e2t


class GeoPySidesCalculator:
    """This algorithm estimates the lengths of the cell sides by calculating
    the positions of the cell faces as the midpoints between centers (using the
    `_extend_from_average` function). It then computes the distance between the
    two opposite faces of each cell.

    Args:
        geodesic (bool): If True, distances are calculated along the
        geodesic (the shortest path between two points on a curved
        surface) on the WGS84 ellipsoid. If False, the Earth is
        approximated as a sphere, and the length of the great circle
        is returned.
    """

    def __init__(self, geodesic: bool = True):
        self._geodesic = bool(geodesic)

    def __call__(self, xlevels, ylevels):
        u_faces_lon_coords = extend_from_average(xlevels, axis=1)
        u_faces_lat_coords = extend_from_average(ylevels, axis=1)
        v_faces_lon_coords = extend_from_average(xlevels, axis=0)
        v_faces_lat_coords = extend_from_average(ylevels, axis=0)

        # For reference, here we report the original Fortran names
        # used inside NEMO model:
        # u_faces_lon_coords = glamu
        # u_faces_lat_coords = gphiu
        # v_faces_lon_coords = glamv
        # v_faces_lat_coords = gphiv

        if self._geodesic:
            distance_function = _geodetic_distance
        else:
            distance_function = _great_circle_distance

        e1t = np.empty_like(xlevels, dtype=np.float32)
        e2t = np.empty_like(xlevels, dtype=np.float32)
        for i, j in cart_prod(range(e1t.shape[0]), range(e1t.shape[1])):
            e1t[i, j] = distance_function(
                (u_faces_lat_coords[i, j], u_faces_lon_coords[i, j]),
                (u_faces_lat_coords[i, j + 1], u_faces_lon_coords[i, j + 1]),
            )
            e2t[i, j] = distance_function(
                (v_faces_lat_coords[i, j], v_faces_lon_coords[i, j]),
                (v_faces_lat_coords[i + 1, j], v_faces_lon_coords[i + 1, j]),
            )

        return e1t, e2t
