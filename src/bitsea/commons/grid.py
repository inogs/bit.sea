from abc import ABC
from abc import abstractmethod
from os import PathLike
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import netCDF4
import numpy as np
from numpy.typing import ArrayLike
from sklearn.neighbors import BallTree

from bitsea.basins.region import Polygon
from bitsea.basins.region import Rectangle
from bitsea.basins.region import Region
from bitsea.commons.geodistances import extend_from_average
from bitsea.commons.geodistances import GeoPySidesCalculator
from bitsea.commons.geodistances import GridSidesAlgorithm
from bitsea.commons.geodistances import NemoGridSidesCalculator
from bitsea.utilities.array_wrapper import BooleanArrayWrapper


class OutsideDomain(ValueError):
    pass


class Regular(ABC):
    """We say that a Grid object is regular if all columns of `xlevels` and all
    rows of `ylevels` are identical. This feature allows the grid's position to
    be represented using only two 1D arrays, resulting in substantial memory
    savings.

    This class extends the interface of any class by providing two
    additional methods for retrieving these 1D arrays.
    """

    @property
    @abstractmethod
    def lon(self):
        """Returns a 1D array with the longitudes of the center points of the
        grid."""
        raise NotImplementedError

    @property
    @abstractmethod
    def lat(self):
        """Returns a 1D array with the latitudes of the center points of the
        grid."""
        raise NotImplementedError


class Grid(ABC):
    """A `Grid` represents a 2D structure that defines the positions of the
    cell centers in our models.

    This class serves as an interface, consolidating all the relevant
    information and functionality that a grid-like object might provide.
    """

    @property
    @abstractmethod
    def xlevels(self) -> np.ndarray:
        """Returns the longitudes of the points of the grid as a 2D numpy
        array."""
        raise NotImplementedError

    @property
    @abstractmethod
    def ylevels(self) -> np.ndarray:
        """Returns the latitudes of the points of the grid as a 2D numpy
        array."""
        raise NotImplementedError

    @property
    @abstractmethod
    def shape(self) -> Tuple[int, ...]:
        """Returns the shape of this object."""
        raise NotImplementedError

    @property
    @abstractmethod
    def coordinate_dtype(self) -> np.dtype:
        """Returns the dtype used to store the coordinates."""
        raise NotImplementedError

    @property
    @abstractmethod
    def e1t(self) -> np.ndarray:
        """A 2D Numpy array such that `e1t[i, j]` is the length of the segment
        of the parallel passing through the cell centered at `(xlevels[i, j],
        ylevels[i, j])` that lies within the cell (in meters)."""
        raise NotImplementedError

    @property
    @abstractmethod
    def e2t(self):
        """A 2D Numpy array such that `e2t[i, j]` is the length of the segment
        of the meridian passing through the cell centered at `(xlevels[i, j],
        ylevels[i, j])` that lies within the cell (in meters)."""
        raise NotImplementedError

    @property
    def area(self):
        """A 2D numpy array such that area[i, j] is the area of the cell `(i,
        j)` (in square meters)."""
        return self.e1t * self.e2t

    @abstractmethod
    def is_inside_domain(
        self,
        *,
        lon: Union[float, ArrayLike],
        lat: Union[float, ArrayLike],
        raise_if_outside: bool = False,
    ):
        """
        Check if one or more points described by their coordinates are over
        one of the cells of the grid.

        Args:
            lon (float or array): Longitude of the point, or an array
              of longitudes if multiple points are being checked.
            lat (float or array): Latitude of the point, or an array
              of latitudes if multiple points are being checked.
            raise_if_outside (bool): If True, raises an `OutsideDomain`
              exception if any point is outside the rectangle.

        Returns:
            An array of boolean

        Raises:
            An `OutsideDomain` exception if the `raise_if_outside` flag is True
            and any point is outside the area occupied by the grid
        """
        raise NotImplementedError

    @abstractmethod
    def convert_lon_lat_to_indices(
        self, *, lon: Union[float, np.ndarray], lat: Union[float, np.ndarray]
    ) -> Tuple:
        """Converts longitude and latitude to the nearest indices on the mask.

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
        Returns:
            a tuple with two numbers i and j, so that `xlevels[i, j]`
            is the longitude of the closest center point of the grid
            and `ylevels[i, j]` is the latitude of the same point.
        """
        raise NotImplementedError

    @abstractmethod
    def convert_i_j_to_lon_lat(
        self, i: Union[int, np.ndarray], j: Union[int, np.ndarray]
    ) -> Tuple:
        """Given two indices, return the corresponding longitude and latitude
        of the point identified by indices i and j.

        Args:
            i (int): The index representing the row.
            j (int): The index representing the column.

        Returns:
            tuple: A tuple containing the longitude and latitude of the point.
        """
        raise NotImplementedError

    def is_regular(self):
        """Returns "True" if the grid is regular, i.e. if all the columns of
        xlevels and all the rows of ylevels are the same."""
        return False

    @staticmethod
    def from_file(
        file_path: PathLike,
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
    ):
        """Reads a NetCDF file and returns a `Grid` object.

        Args:
            file_path (PathLike): Path to the NetCDF file.
            ylevels_var_name (str): Name of the variable in the NetCDF
              file that contains the y-levels (latitude values).
            xlevels_var_name (str): Name of the variable in the NetCDF
              file that contains the x-levels (longitude values).

        Returns:
            Grid: A `Grid` object. If the file contains a regular grid,
            a `RegularGrid` object is returned; otherwise it produces
            an IrregularGrid.
        """
        with netCDF4.Dataset(file_path) as f:
            return Grid.from_file_pointer(f, ylevels_var_name, xlevels_var_name)

    @staticmethod
    def from_file_pointer(
        file_pointer: netCDF4.Dataset,
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
    ):
        """Reads a NetCDF file and returns a `Grid` object.

        This function is very similar to the `from_file` static method;
        the only difference is that this function requires a file_pointer
        instead of the Path of the file

        Args:
            file_pointer (netCDF4.Dataset): A open NetCDF Dataset.
            ylevels_var_name (str): Name of the variable in the NetCDF
              file that contains the y-levels (latitude values).
            xlevels_var_name (str): Name of the variable in the NetCDF
              file that contains the x-levels (longitude values).

        Returns:
            Grid: A `Grid` object. If the file contains a regular grid,
            a `RegularGrid` object is returned; otherwise it produces
            an IrregularGrid.
        """

        ylevels = np.asarray(file_pointer.variables[ylevels_var_name])
        xlevels = np.asarray(file_pointer.variables[xlevels_var_name])

        if len(xlevels.shape) == 4:
            xlevels = xlevels[0, 0, :, :]
        if len(ylevels.shape) == 4:
            ylevels = ylevels[0, 0, :, :]

        if "e1t" in file_pointer.variables:
            e1t = np.asarray(
                file_pointer.variables["e1t"][0, 0, :, :], dtype=np.float32
            )
        else:
            e1t = None

        if "e2t" in file_pointer.variables:
            e2t = np.asarray(
                file_pointer.variables["e2t"][0, 0, :, :], dtype=np.float32
            )
        else:
            e2t = None

        # Now we check if the grid is regular by checking if replicating
        # the first row of the longitudes and the first column of the
        # latitudes we get the same grid that we have read from the file
        x1d = xlevels[0, :]
        y1d = ylevels[:, 0]
        x2d = x1d[np.newaxis, :]
        y2d = y1d[:, np.newaxis]
        dist = np.max((x2d - xlevels) ** 2 + (y2d - ylevels) ** 2)
        if dist < 1e-8:
            return RegularGrid(lon=x1d, lat=y1d, e1t=e1t, e2t=e2t)

        return IrregularGrid(xlevels=xlevels, ylevels=ylevels, e1t=e1t, e2t=e2t)


class BaseGrid(Grid, ABC):
    """
    The `BaseGrid` class represents a generic grid structure capable of
    calculating vector dimensions that represent the grid size in meters.

    This base class provides a shared algorithm for computing grid sizes,
    used by both `RegularGrid` and `IrregularGrid`, making it a common
    ancestor for these grid types.
    """

    def __init__(
        self,
        *,
        e1t: Optional[np.ndarray] = None,
        e2t: Optional[np.ndarray] = None,
        **kwargs,
    ):
        if e1t is not None:
            self._e1t = e1t.view()
            self._e1t.setflags(write=False)
        else:
            self._e1t = None

        if e2t is not None:
            self._e2t = e2t.view()
            self._e2t.setflags(write=False)
        else:
            self._e2t = None

        self._default_side_algorithm = GridSidesAlgorithm.GEODESIC

        self._domain: Optional[Region] = None

        # Mixin behaviour
        super().__init__(**kwargs)

    @abstractmethod
    def _build_domain(self) -> Region:
        raise NotImplementedError

    def is_inside_domain(
        self,
        *,
        lon: Union[float, ArrayLike],
        lat: Union[float, ArrayLike],
        raise_if_outside: bool = False,
    ):
        if self._domain is None:
            self._domain = self._build_domain()
        output = self._domain.is_inside(lon=lon, lat=lat)

        if not raise_if_outside or np.all(output):
            return output

        if output.ndim > 0:
            outside_indices = np.nonzero(~output)
            first_point = tuple(p[0] for p in outside_indices)
            raise OutsideDomain(
                f"Point {first_point} with coords lon = "
                f"{lon[first_point]} and lat = {lat[first_point]} "
                f"is outside the domain"
            )

        raise OutsideDomain(
            f"Point with coords (lon = {lon} and lat = {lat}) is "
            f"outside the domain"
        )

    def _build_e_t_arrays(
        self, algorithm: GridSidesAlgorithm
    ) -> Tuple[np.ndarray, np.ndarray]:
        xlevels = self.xlevels
        ylevels = self.ylevels

        if algorithm == GridSidesAlgorithm.NEMO:
            calculator = NemoGridSidesCalculator()
        elif algorithm == GridSidesAlgorithm.GREAT_CIRCLE:
            calculator = GeoPySidesCalculator(geodesic=False)
        elif algorithm == GridSidesAlgorithm.GEODESIC:
            calculator = GeoPySidesCalculator(geodesic=True)
        else:
            raise ValueError("Not implemented algorithm: {}".format(algorithm))

        return calculator(xlevels=xlevels, ylevels=ylevels)

    def _build_and_store_e_t_arrays(self):
        e1t, e2t = self._build_e_t_arrays(
            algorithm=self._default_side_algorithm
        )
        if self._e1t is None:
            self._e1t = e1t
            self._e1t.setflags(write=False)
        if self._e2t is None:
            self._e2t = e2t
            self._e2t.setflags(write=False)


class IrregularGrid(BaseGrid):
    """Represents a 2D discretization of a geographical region, storing the
    positions of cell centers in terms of longitude and latitude.

    Args:
        xlevels (2D array): Longitudes of the cell centers.
        ylevels (2D array): Latitudes of the cell centers.
        e1t (2D array): e1t[i, j] must be the length of the segment of
          the parallel passing through the cell `(i, j)` that lies
          within the cell. If this information is not available, this
          object can be `None`, and an approximation will be computed
          by this object.
        e2t (2D array): is like e1t but along the meridian. It also can
          be `None`.
    """

    def __init__(
        self,
        *,
        xlevels: np.ndarray,
        ylevels: np.ndarray,
        e1t: Optional[np.ndarray] = None,
        e2t: Optional[np.ndarray] = None,
    ):
        xlevels = np.asarray(xlevels)
        ylevels = np.asarray(ylevels)

        try:
            xlevels, ylevels = np.broadcast_arrays(xlevels, ylevels)
        except ValueError:
            raise ValueError(
                f"The longitude and latitude arrays are expected to have the "
                f"same shape (and both of them should be 2D arrays): "
                f"currently lon = {xlevels.shape} and lat = {ylevels.shape}"
            )

        if len(xlevels.shape) != 2:
            raise ValueError(
                f"The lon array is expected to be a 2D array; its current "
                f"shape is {xlevels.shape}"
            )

        if xlevels.dtype != ylevels.dtype:
            raise ValueError(
                f"The longitude and latitude arrays are expected to have the "
                f"same dtype: currently lon = {xlevels.dtype} and lat = "
                f"{ylevels.dtype}"
            )

        self._xlevels = xlevels.view()
        self._ylevels = ylevels.view()

        self._xlevels.setflags(write=False)
        self._ylevels.setflags(write=False)

        super().__init__(e1t=e1t, e2t=e2t)

        # BallTree to find the nearest cell center to a given point.
        # Initialized as a list that contains only `None`; when we need a
        # BallTree we will replace the first element of the list with a BallTree
        # object (the list will never have more than one element).
        # This mechanism allows to share the cache among different objects that
        # represent the same grid (see, for example, the MaskLayers objects);
        # when one of the objects initialize the cache, it is initialized for
        # all of them
        self._balltree_cache: List[Optional[BallTree]] = [None]

    @property
    def xlevels(self) -> np.ndarray:
        return self._xlevels

    @property
    def ylevels(self) -> np.ndarray:
        return self._ylevels

    @property
    def shape(self) -> Tuple[int, ...]:
        return self._xlevels.shape

    @property
    def coordinate_dtype(self):
        return self._xlevels.dtype

    @property
    def lon(self):
        raise ValueError(
            "This grid is not regular, so the 1D lon array is not defined. "
            "Use .xlevels instead."
        )

    @property
    def lat(self):
        raise ValueError(
            "This grid is not regular, so the 1D lat array is not defined. "
            "Use .ylevels instead."
        )

    @property
    def e1t(self) -> np.ndarray:
        if self._e1t is None:
            self._build_and_store_e_t_arrays()
        return self._e1t

    @property
    def e2t(self) -> np.ndarray:
        if self._e2t is None:
            self._build_and_store_e_t_arrays()
        return self._e2t

    def _initialize_balltree(self):
        """Initializes the ball tree for this object, enabling efficient
        nearest-neighbor searches."""
        # Here we put ylevels and then xlevels because the balltree
        # needs the first coordinate to be the latitude
        data = np.stack((self.ylevels, self.xlevels), axis=-1)
        data = np.radians(data.reshape(-1, 2))

        # noinspection PyArgumentList
        self._balltree_cache[0] = BallTree(data, metric="haversine")

    def _build_domain(self):
        # We initialize a Polygon by attaching together all the center of
        # the faces that are on the boundary of the domain. Here we compute
        # the coordinates of the faces
        u_faces_lon_coords = extend_from_average(self.xlevels, axis=1)
        u_faces_lat_coords = extend_from_average(self.ylevels, axis=1)
        v_faces_lon_coords = extend_from_average(self.xlevels, axis=0)
        v_faces_lat_coords = extend_from_average(self.ylevels, axis=0)

        # The points of the boundary of the polygon that are on the right
        # The first and last vertices are choosen accordingly to the lower and
        # upper value of the faces on the top and bottom
        right_boundary_lon = np.concatenate(
            (
                (u_faces_lon_coords[0, -1],),
                u_faces_lon_coords[:, -1],
                (u_faces_lon_coords[-1, -1],),
            )
        )
        right_boundary_lat = np.concatenate(
            (
                (v_faces_lat_coords[0, -1],),
                u_faces_lat_coords[:, -1],
                (v_faces_lat_coords[-1, -1],),
            )
        )

        # Now we do the same for the other sides. Here we need to concatenate
        # only two elements because one of the vertex is shared with the
        # previous side
        top_boundary_lon = np.concatenate(
            (v_faces_lon_coords[-1, :][::-1], (u_faces_lon_coords[-1, 0],))
        )
        top_boundary_lat = np.concatenate(
            (v_faces_lat_coords[-1, :][::-1], (v_faces_lat_coords[-1, 0],))
        )

        left_boundary_lon = np.concatenate(
            (u_faces_lon_coords[:, 0][::-1], (u_faces_lon_coords[0, 0],))
        )
        left_boundary_lat = np.concatenate(
            (u_faces_lat_coords[:, 0][::-1], (v_faces_lat_coords[0, 0],))
        )

        bottom_boundary_lon = v_faces_lon_coords[0, :]
        bottom_boundary_lat = v_faces_lat_coords[0, :]

        lon_boundary = np.concatenate(
            (
                right_boundary_lon,
                top_boundary_lon,
                left_boundary_lon,
                bottom_boundary_lon,
            )
        )
        lat_boundary = np.concatenate(
            (
                right_boundary_lat,
                top_boundary_lat,
                left_boundary_lat,
                bottom_boundary_lat,
            )
        )

        return Polygon(lon_list=lon_boundary, lat_list=lat_boundary)

    def convert_lon_lat_to_indices(
        self, *, lon: Union[float, np.ndarray], lat: Union[float, np.ndarray]
    ) -> Tuple:
        # If lon and lat are numpy array, check if their shape is compatible
        # by simply trying a broadcast
        np.broadcast(lon, lat)

        # Raise an OutsideDomain error if a point is too far from the grid
        self.is_inside_domain(lon=lon, lat=lat, raise_if_outside=True)

        lon, lat = np.broadcast_arrays(lon, lat)

        # Ensure that the balltree is initialized
        if self._balltree_cache[0] is None:
            self._initialize_balltree()
        balltree = self._balltree_cache[0]

        # Reshape the points so that it becomes a couple of longitudes and
        # latitudes
        query_data = np.radians(np.stack((lat, lon), axis=-1))
        assert query_data.shape[-1] == 2

        # BallTree wants two-dimensional arrays
        if len(query_data.shape) == 1:
            query_data = query_data[np.newaxis, :]
        if len(query_data.shape) > 2:
            query_data = query_data.reshape((-1, 2))

        _, center_indices = balltree.query(query_data)
        center_indices = center_indices.reshape(lon.shape)

        # We need to reshape again the indices to be compatible with the
        # shape of xlevels and ylevels
        i_indices = center_indices // self.xlevels.shape[-1]
        j_indices = center_indices % self.xlevels.shape[-1]

        if not isinstance(i_indices, np.ndarray):
            i_indices = int(i_indices)
        if not isinstance(j_indices, np.ndarray):
            j_indices = int(j_indices)

        return i_indices, j_indices

    def convert_i_j_to_lon_lat(
        self, i: Union[int, np.ndarray], j: Union[int, np.ndarray]
    ) -> Tuple:
        return self._xlevels[i, j], self._ylevels[i, j]


class RegularGrid(BaseGrid, Regular):
    """
    A `RegularGrid` is defined by two vectors: `lon` (longitude) and
    `lat` (latitude).
    It behaves like a standard `Grid`, where the grid points are
    defined as:

        - `xlevels[i, j] = lon[j]` for every `i`
        - `ylevels[i, j] = lat[i]` for every `j`

    Args:
        lon (1D array): Longitudes of the cell centers.
        lat (1D array): Latitudes of the cell centers.
        e1t (2D array): e1t[i, j] must be the length of the segment of
          the parallel passing through the cell `(i, j)` that lies
          within the cell. If this information is not available, this
          object can be `None`, and an approximation will be computed
          by this object.
        e2t (2D array): is like e1t but along the meridian. It also can
          be `None`.
    """

    def __init__(
        self,
        *,
        lon,
        lat,
        e1t: Optional[np.ndarray] = None,
        e2t: Optional[np.ndarray] = None,
    ):
        lon = np.asarray(lon).view()
        lat = np.asarray(lat).view()

        if len(lon.shape) != 1:
            raise ValueError(
                f"lon should be a 1D array: its current shape is {lon.shape}"
            )
        if len(lat.shape) != 1:
            raise ValueError(
                f"lat should be a 1D array: its current shape is {lat.shape}"
            )

        self._lon = lon
        self._lat = lat

        self.lon.setflags(write=False)
        self.lat.setflags(write=False)

        if np.any(self._lon[1:] - self._lon[:-1] <= 0):
            raise ValueError(
                f"lon vector must be strictly increasing; received {self._lon}"
            )
        if np.any(self._lat[1:] - self._lat[:-1] <= 0):
            raise ValueError(
                f"lat vector must be strictly increasing; received {self._lon}"
            )

        if self._lon.dtype != self._lat.dtype:
            raise ValueError(
                f"The longitude and latitude arrays are expected to have the "
                f"same dtype: currently lon = {self._lon.dtype} and lat = "
                f"{self._lat.dtype}"
            )

        super().__init__(e1t=e1t, e2t=e2t)
        self._default_side_algorithm = GridSidesAlgorithm.NEMO

    @property
    def lon(self):
        return self._lon

    @property
    def lat(self):
        return self._lat

    @property
    def xlevels(self):
        grid_shape = (self._lat.shape[0], self._lon.shape[0])
        return np.broadcast_to(self._lon, grid_shape)

    @property
    def ylevels(self):
        grid_shape = (self._lat.shape[0], self._lon.shape[0])
        return np.broadcast_to(self._lat[:, np.newaxis], grid_shape)

    @property
    def shape(self) -> Tuple[int, ...]:
        return self._lat.shape[0], self._lon.shape[0]

    @property
    def coordinate_dtype(self):
        return self._lon.dtype

    @property
    def e1t(self) -> np.ndarray:
        if self._e1t is None:
            self._build_and_store_e_t_arrays()
        return self._e1t

    @property
    def e2t(self) -> np.ndarray:
        if self._e2t is None:
            self._build_and_store_e_t_arrays()
        return self._e2t

    def is_regular(self):
        return True

    def convert_i_j_to_lon_lat(self, i, j):
        return self.lon[j], self.lat[i]

    def _build_domain(self):
        faces_lon = extend_from_average(self.lon)
        faces_lat = extend_from_average(self.lat)
        return Rectangle(
            lonmin=faces_lon[0],
            lonmax=faces_lon[-1],
            latmin=faces_lat[0],
            latmax=faces_lat[-1],
        )

    def convert_lon_lat_to_indices(self, *, lon, lat):
        return_array = False
        if isinstance(lon, np.ndarray) or isinstance(lat, np.ndarray):
            return_array = True

        lon, lat = np.broadcast_arrays(lon, lat)

        # Raise an OutsideDomain error if a point is too far from the grid
        self.is_inside_domain(lon=lon, lat=lat, raise_if_outside=True)

        lon_distances = np.abs(self.lon[:, np.newaxis] - lon)
        lat_distances = np.abs(self.lat[:, np.newaxis] - lat)

        i = np.argmin(lon_distances, axis=0)
        j = np.argmin(lat_distances, axis=0)

        if not return_array:
            i = np.squeeze(i)
            j = np.squeeze(j)
            assert i.shape == tuple()
            assert j.shape == tuple()
            i = int(i.item())
            j = int(j.item())

        return i, j

    @staticmethod
    def from_file(
        file_name: PathLike,
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
    ):
        grid = IrregularGrid.from_file(
            file_name,
            ylevels_var_name=ylevels_var_name,
            xlevels_var_name=xlevels_var_name,
        )

        if not grid.is_regular():
            raise ValueError(
                f"File {file_name} does not contain a regular grid"
            )

        return grid


class MaskLayer(BooleanArrayWrapper, Grid, ABC):
    """The `MaskLayer` class represents a single layer of a mask at a defined
    depth. It encapsulates all information from a `Grid`, along with additional
    attributes: the depth, the thickness of the layer, and a boolean mask.
    The mask indicates whether each cell contains water (`True`) or lies on
    land (`False`)."""

    def __init__(
        self,
        *,
        depth: float,
        thickness: float,
        mask: np.ndarray,
        e1t: Optional[np.ndarray] = None,
        e2t: Optional[np.ndarray] = None,
        **kwargs,
    ):
        self._depth = depth
        self._thickness = thickness

        super().__init__(wrapped_data=mask, e1t=e1t, e2t=e2t, **kwargs)

        if (
            super(BooleanArrayWrapper.__bases__[0], self).shape
            != self._data_array.shape
        ):
            raise ValueError(
                "The mask array is not coherent with the current grid"
            )

    @property
    def depth(self) -> float:
        return self._depth

    @property
    def thickness(self) -> float:
        return self._thickness

    @classmethod
    def from_grid(
        cls, grid: Grid, *, depth: float, thickness: float, mask: np.ndarray
    ):
        """
        Extends a grid object and transforms it into a grid. This is usually
        more
        efficient than using the initializer directly, as it can reuse some cached
        values from the original grid object, reducing the need for copying.
        """
        if isinstance(grid, RegularGrid):
            return RegularMaskLayer.from_grid(
                grid, depth=depth, thickness=thickness, mask=mask
            )
        else:
            return IrregularMaskLayer.from_grid(
                grid, depth=depth, thickness=thickness, mask=mask
            )


class RegularMaskLayer(MaskLayer, RegularGrid):
    def __init__(
        self,
        *,
        lon: np.ndarray,
        lat: np.ndarray,
        depth: float,
        thickness: float,
        mask: np.ndarray,
        e1t: Optional[np.ndarray] = None,
        e2t: Optional[np.ndarray] = None,
    ):
        super().__init__(
            depth=depth,
            thickness=thickness,
            mask=mask,
            lon=lon,
            lat=lat,
            e1t=e1t,
            e2t=e2t,
        )

    @classmethod
    def from_grid(
        cls, grid: Grid, *, depth: float, thickness: float, mask: np.ndarray
    ):
        if not isinstance(grid, RegularGrid):
            raise ValueError("grid must be a RegularGrid")

        return cls(
            lon=grid.lon,
            lat=grid.lat,
            depth=depth,
            thickness=thickness,
            mask=mask,
            e1t=grid._e1t,
            e2t=grid._e2t,
        )


class IrregularMaskLayer(MaskLayer, IrregularGrid):
    def __init__(
        self,
        *,
        xlevels: np.ndarray,
        ylevels: np.ndarray,
        depth: float,
        thickness: float,
        mask: np.ndarray,
        e1t: Optional[np.ndarray] = None,
        e2t: Optional[np.ndarray] = None,
    ):
        super().__init__(
            depth=depth,
            thickness=thickness,
            mask=mask,
            xlevels=xlevels,
            ylevels=ylevels,
            e1t=e1t,
            e2t=e2t,
        )

    @classmethod
    def from_grid(
        cls, grid: Grid, *, depth: float, thickness: float, mask: np.ndarray
    ):
        if not isinstance(grid, IrregularGrid):
            raise ValueError("grid must be a IrregularGrid")

        output = cls(
            xlevels=grid.xlevels,
            ylevels=grid.ylevels,
            depth=depth,
            thickness=thickness,
            mask=mask,
            e1t=grid._e1t,
            e2t=grid._e2t,
        )

        output._balltree_cache = grid._balltree_cache

        return output
