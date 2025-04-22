from numbers import Real
from os import PathLike
from typing import List
from typing import Literal
from typing import Optional
from typing import Tuple
from typing import Union
from warnings import warn

import netCDF4
import numpy as np
from numpy.typing import ArrayLike

from bitsea.commons.geodistances import extend_from_average
from bitsea.commons.grid import Grid
from bitsea.commons.grid import IrregularGrid
from bitsea.commons.grid import RegularGrid


class Mesh:
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
    def jpk(self) -> int:
        return self._zlevels.shape[0]

    @property
    def jpj(self) -> int:
        return self._grid.shape[0]

    @property
    def jpi(self) -> int:
        return self._grid.shape[1]

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
    def area(self):
        """A 2D numpy array such that area[i, j] is the area of the cell `(i,
        j)` (in square meters)."""
        return self.e1t * self.e2t

    @property
    def lon(self):
        try:
            # noinspection PyUnresolvedReferences
            return self._grid.lon
        except AttributeError:
            raise AttributeError(
                f"This {self.__class__.__name__} object is constructed on a "
                f'Grid that is not a RegularGrid; therefore, the "lon" '
                f"attribute is not available"
            )

    @property
    def lat(self):
        try:
            # noinspection PyUnresolvedReferences
            return self._grid.lat
        except AttributeError:
            raise AttributeError(
                f"This {self.__class__.__name__} object is constructed on a "
                f'Grid that is not a RegularGrid; therefore, the "lon" '
                f"attribute is not available"
            )

    def is_inside_domain(
        self,
        *,
        lon: Union[float, ArrayLike],
        lat: Union[float, ArrayLike],
        raise_if_outside: bool = False,
    ):
        return self._grid.is_inside_domain(
            lon=lon, lat=lat, raise_if_outside=raise_if_outside
        )

    @property
    def dz(self) -> np.ndarray:
        """`dz` is a 1D version of the `e3t` array, representing the vertical
        size of the cells.

        It does not account for the seabed, so each depth level has a
        single value, making it a 1D array.
        """
        return np.max(self.e3t, axis=(1, 2))

    def convert_lat_lon_to_indices(
        self, *, lon: Union[float, np.ndarray], lat: Union[float, np.ndarray]
    ) -> Tuple:
        return self._grid.convert_lat_lon_to_indices(lon=lon, lat=lat)

    def convert_lon_lat_to_indices(
        self, *, lon: Union[float, np.ndarray], lat: Union[float, np.ndarray]
    ) -> Tuple:
        return self._grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    def convert_i_j_to_lon_lat(
        self, i: Union[int, np.ndarray], j: Union[int, np.ndarray]
    ) -> Tuple:
        return self._grid.convert_i_j_to_lon_lat(i=i, j=j)

    def getDepthIndex(self, z: Real) -> int:
        warn(
            'This method is deprecated. Use "get_depth_index" instead.',
            DeprecationWarning,
        )
        return self.get_depth_index(z)

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

    def column_side_area(
        self,
        ji: int,
        jj: int,
        side: Literal["N", "S", "W", "E"],
        n_vertical_cells: int,
    ) -> float:
        """
        Calculates the lateral area of a water column with the specified depth.

        Args:
            ji (int): Longitudinal index.
            jj (int): Latitudinal index.
            side (Literal["N", "S", "W", "E"]): The side of the column for
              which to compute the area, specified by the first letter of a
              cardinal direction.
            n_vertical_cells (int): The number of cells in the column.

        Returns:
            float: The lateral area of the specified column side.
        """
        if n_vertical_cells is None:
            raise ValueError(
                "A `mesh` object does not know which cells contain water. Therefore, "
                "the `column_side_area` method requires an explicit number of vertical "
                "cells; the parameter `n_vertical_cells` can not be `None`"
            )

        if side in ["E", "W"]:
            return self.e2t[jj, ji] * self.e3t[:n_vertical_cells, jj, ji].sum()
        elif side in ["N", "S"]:
            return self.e1t[jj, ji] * self.e3t[:n_vertical_cells, jj, ji].sum()
        else:
            raise ValueError(f'Invalid side: "{side}"')

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
        *,
        zlevels_var_name: str = "nav_lev",
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
        e3t_var_name: Optional[str] = None,
        read_e3t: bool = False,
    ):
        with netCDF4.Dataset(file_path, "r") as f:
            return cls.from_file_pointer(
                f,
                zlevels_var_name=zlevels_var_name,
                ylevels_var_name=ylevels_var_name,
                xlevels_var_name=xlevels_var_name,
                e3t_var_name=e3t_var_name,
                read_e3t=read_e3t,
            )

    @classmethod
    def from_file_pointer(
        cls,
        file_pointer: netCDF4.Dataset,
        *,
        zlevels_var_name: str = "nav_lev",
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
        e3t_var_name: Optional[str] = None,
        read_e3t: bool = False,
    ):
        grid = Grid.from_file_pointer(
            file_pointer,
            ylevels_var_name=ylevels_var_name,
            xlevels_var_name=xlevels_var_name,
        )

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
            if e3t_var_name is None:
                e3t_names = ("e3t", "e3t_0")
            else:
                e3t_names = (e3t_var_name,)
            for e3t_name in e3t_names:
                if e3t_name in file_pointer.variables:
                    e3t = np.ma.getdata(
                        file_pointer.variables[e3t_name][0, :, :, :]
                    )
                    break
            else:
                raise KeyError(
                    "No e3t variable found inside netCDF file; "
                    f"tried the following names: {e3t_names}"
                )
        else:
            e3t = None

        return cls(grid, zlevels, e3t)

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
        return Mesh(grid, zlevels)
