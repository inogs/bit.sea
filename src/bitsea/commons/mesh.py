from numbers import Real
from os import PathLike
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import netCDF4
import numpy as np
from numpy.typing import ArrayLike

from bitsea.commons.geodistances import extend_from_average
from bitsea.commons.grid import Grid
from bitsea.commons.grid import IrregularGrid
from bitsea.commons.grid import Regular
from bitsea.commons.grid import RegularGrid


class Mesh(Grid):
    """A `Mesh` is a 3D extension of a `Grid`, containing all its information,
    while also adding details about vertical levels.

    Args:
        grid (Grid): A `Grid` providing 2D grid information.
        zlevels (np.ndarray): A 1D array representing the depth (in
          meters) of cell centers along the vertical axis.
        e3t (Optional[np.ndarray]): A 3D array representing the
          vertical size (in meters) of each cell. This accounts for
          variations in cell size, especially at the bottom layer where
          cells may reach the seabed.

          If `None`, the vertical sizes (`e3t`) are calculated by assuming
          `zlevels` are positioned in the center of each cell. In this
          case, `e3t[:, i, j]` will be uniform across the grid, and
          detailed information about the bottom layer may be lost.
    """

    def __init__(
        self, grid: Grid, zlevels: ArrayLike, e3t: Optional[np.ndarray] = None
    ):
        self._grid = grid

        self._zlevels = np.asarray(zlevels).view()
        self._zlevels.setflags(write=False)

        if self._zlevels.ndim != 1:
            raise ValueError(
                f"zlevels must be a 1D array: its current shape is "
                f"{self._zlevels.shape}"
            )

        if np.any(self._zlevels < 0):
            negative_index = np.where(self._zlevels < 0)[0]
            raise ValueError(
                f"zlevels must be non-negative, but zlevels[{negative_index}] "
                f"is {zlevels[negative_index]}"
            )

        if np.any(self._zlevels[:-1] - self._zlevels[1:] >= 0):
            raise ValueError("zlevels values should be strictly increasing")

        if e3t is not None:
            expected_shape = (self._zlevels.shape[0],) + self._grid.shape
            e3t = np.asarray(e3t).view()

            self._e3t = e3t
            # If it is a 1 dim array, we try to broadcast it on all the grid
            if self._e3t.ndim == 1:
                self._e3t = self._e3t[:, np.newaxis, np.newaxis]

            if self._e3t.shape != expected_shape:
                try:
                    self._e3t = np.broadcast_to(self._e3t, expected_shape)
                except ValueError:
                    raise ValueError(
                        f"e3t is expected to be {expected_shape}, but got "
                        f"an array with shape {e3t.shape}"
                    )
            self._e3t.setflags(write=False)
        else:
            # We construct e3t so that zlevels is its average
            layer_boundaries = extend_from_average(self._zlevels, 0, 0.0)
            e3t_md = layer_boundaries[1:] - layer_boundaries[:-1]
            self._e3t = np.broadcast_to(
                e3t_md[:, np.newaxis, np.newaxis],
                (self._zlevels.shape[0],) + self._grid.shape,
            )

        super().__init__()

    @property
    def grid(self) -> Grid:
        """Returns a `Grid` that contains the 2D information of this object,
        but not the zlevels."""
        return self._grid

    def is_regular(self):
        return self._grid.is_regular()

    @property
    def xlevels(self) -> np.ndarray:
        return self._grid.xlevels

    @property
    def ylevels(self) -> np.ndarray:
        return self._grid.ylevels

    @property
    def zlevels(self) -> np.ndarray:
        return self._zlevels

    @property
    def shape(self) -> Tuple[int, ...]:
        return (self._zlevels.shape[0],) + self._grid.shape

    @property
    def coordinate_dtype(self) -> np.dtype:
        return self._grid.coordinate_dtype

    @property
    def e1t(self) -> np.ndarray:
        return self._grid.e1t

    @property
    def e2t(self):
        return self._grid.e2t

    @property
    def e3t(self) -> np.ndarray:
        return self._e3t

    @property
    def dz(self) -> np.ndarray:
        """`dz` is a 1D version of the `e3t` array, representing the vertical
        size of the cells.

        It does not account for the seabed, so each depth level has a
        single value, making it a 1D array.
        """
        return np.max(self.e3t, axis=(1, 2))

    def convert_lon_lat_to_indices(
        self, *, lon: Union[float, np.ndarray], lat: Union[float, np.ndarray]
    ) -> Tuple:
        return self._grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    def convert_i_j_to_lon_lat(
        self, i: Union[int, np.ndarray], j: Union[int, np.ndarray]
    ) -> Tuple:
        return self._grid.convert_i_j_to_lon_lat(i=i, j=j)

    def get_depth_index(self, z: Real) -> int:
        """Converts a depth expressed in meters to the corresponding index
        level.

        The returned value is an integer indicating the previous
        (*not* the nearest) depth in the z-levels.
        If `z` is above the first level, 0 is returned.

        Example:
        >>> grid = RegularGrid(np.linspace(0, 10, 11), np.linspace(0, 10, 23))
        >>> mesh = Mesh(grid, zlevels=5 + np.arange(30) * 10.)
        >>> k = mesh.get_depth_index(200.)
        >>> mesh.zlevels[k]
        195.
        """
        z = np.asarray(z)
        output = np.maximum(
            np.searchsorted(self._zlevels, z, side="right") - 1, 0
        )
        if output.ndim == 0:
            output = int(output)
        return output

    @staticmethod
    def from_levels(
        xlevels: np.ndarray, ylevels: np.ndarray, zlevels: np.ndarray
    ):
        """Create a `Mesh` from the three arrays that define the position of
        the centers of the cells. This is the minimum amount of information
        needed to generate a Mesh (the size of the cells in each direction will
        be approximate by the algorithms of `bit.sea`)

        Args:
            xlevels (np.ndarray): A 2D array containing the longitudes
              of the grid cell centers.
            ylevels (np.ndarray): A 2D array containing the latitudes
              of the grid cell centers.
            zlevels (np.ndarray): A 1D array containing the depth (in
              meters) of the centers of the cells.
        """
        grid = IrregularGrid(xlevels=xlevels, ylevels=ylevels)
        return Mesh(grid, zlevels)

    @classmethod
    def from_file(
        cls,
        file_path: PathLike,
        zlevels_var_name: str = "nav_lev",
        read_e3t: bool = False,
    ):
        with netCDF4.Dataset(file_path, "r") as f:
            return cls.from_file_pointer(f, zlevels_var_name, read_e3t)

    @classmethod
    def from_file_pointer(
        cls,
        file_pointer: netCDF4.Dataset,
        zlevels_var_name: str = "nav_lev",
        read_e3t: bool = False,
    ):
        grid = Grid.from_file_pointer(file_pointer)

        zlevels = np.array(
            file_pointer.variables[zlevels_var_name][:], dtype=np.float32
        )

        if zlevels.ndim != 1:
            # Reduce zlevels to a 1D array
            if zlevels.ndim == 2:
                raise ValueError("zlevels can not be a 2D array")

            if zlevels.ndim > 4:
                raise ValueError("zlevels can not have more than 4 dimensions")

            zlevels_slice: List[Union[int, slice]] = [
                0 for _ in range(zlevels.ndim)
            ]
            zlevels_slice[-3] = slice(None)

            zlevels = zlevels[tuple(zlevels_slice)]

        if read_e3t:
            e3t_names = ("e3t", "e3t_0")
            for e3t_name in e3t_names:
                if e3t_name in file_pointer.variables:
                    e3t = np.ma.getdata(
                        file_pointer.variables[e3t_name][0, :, :, :]
                    )
                    break
            else:
                raise KeyError("No e3t variable found inside netCDF file.")
        else:
            e3t = None

        if isinstance(grid, RegularGrid):
            return RegularMesh(grid, zlevels, e3t)

        return cls(grid, zlevels, e3t)


class RegularMesh(Mesh, Regular):
    def __init__(
        self,
        grid: RegularGrid,
        zlevels: ArrayLike,
        e3t: Optional[np.ndarray] = None,
    ):
        if not grid.is_regular():
            raise ValueError("regular_grid argument must be a regular grid.")
        super().__init__(grid, zlevels, e3t)
        self._grid: RegularGrid

    @property
    def lon(self):
        # noinspection PyUnresolvedReferences
        return self._grid.lon

    @property
    def lat(self):
        # noinspection PyUnresolvedReferences
        return self._grid.lat

    @staticmethod
    def from_coordinates(
        *, lon: np.ndarray, lat: np.ndarray, zlevels: np.ndarray
    ):
        """Create a `Mesh` from the three arrays that define the position of
        the centers of the cells. This function is very similar to the
        `from_levels` method of the class `Mesh`, but the coordinates of the
        regular grid are defined by two 1D arrays.

        Args:
            lon (np.ndarray): A 1D array containing the longitudes
              of the grid cell centers.
            lat (np.ndarray): A 1D array containing the latitudes
              of the grid cell centers.
            zlevels (np.ndarray): A 1D array containing the depth (in
              meters) of the centers of the cells.
        """
        grid = RegularGrid(lon=lon, lat=lat)
        return RegularMesh(grid, zlevels)
