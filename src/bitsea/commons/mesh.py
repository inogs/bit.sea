from typing import Optional
from typing import Tuple
from typing import Union

import numpy as np
from numpy.typing import ArrayLike

from bitsea.commons.grid import GridDescriptor, RegularGridDescriptor
from bitsea.commons.grid import extend_from_average


class Mesh(GridDescriptor):
    """
    A `Mesh` is a 3d extension of a `Grid`.
    """
    def __init__(self, grid: GridDescriptor, zlevels: np.ndarray,
                 e3t: Optional[np.ndarray] = None):
        self._grid = grid

        self._zlevels = np.asarray(zlevels).view()
        self._zlevels.setflags(write=False)

        if len(self._zlevels.shape) != 1:
            raise ValueError(
                f'zlevels must be a 1D array: its current shape is '
                f'{self._zlevels.shape}'
            )

        if np.any(self._zlevels < 0):
            negative_index = np.where(self._zlevels < 0)[0]
            raise ValueError(
                f'zlevels must be non-negative, but zlevels[{negative_index}] '
                f'is {zlevels[negative_index]}'
            )

        if np.any(self._zlevels[:-1] - self._zlevels[1:] >= 0):
            raise ValueError(
                'zlevels values should be strictly increasing'
            )

        if e3t is not None:
            expected_shape = (self._zlevels.shape[0],) + self._grid.shape
            self._e3t = np.asarray(e3t).view()
            if self._e3t.shape != expected_shape:
                raise ValueError(
                    f'e3t is expected to be {expected_shape}, but got '
                    f'an array with shape {self._e3t.shape}'
                )
            self._e3t.setflags(write=False)
        else:
            layer_boundaries = extend_from_average(self._zlevels, 0, 0.)
            e3t_md = (layer_boundaries[1:] + layer_boundaries[:-1]) / 2.
            self._e3t = np.broadcast_to(
                e3t_md[:, np.newaxis, np.newaxis],
                (self._zlevels.shape[0],) + self._grid.shape
            )

    @property
    def grid(self):
        return self._grid

    @property
    def is_regular(self):
        return self._grid.is_regular

    @property
    def xlevels(self) -> np.ndarray:
        return self._grid.xlevels

    @property
    def ylevels(self) -> np.ndarray:
        return self._grid.ylevels

    def zlevels(self):
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
    def dz(self):
        return np.max(self.e3t, axis=(1, 2))

    def convert_lon_lat_to_indices(self, *, lon: Union[float, np.ndarray],
                                   lat: Union[float, np.ndarray]) -> Tuple:
        return self._grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    def convert_i_j_to_lon_lat(self, i: Union[int, np.ndarray],
                               j: Union[int, np.ndarray]) -> Tuple:
        return self._grid.convert_i_j_to_lon_lat(i=i, j=j)

    def get_depth_index(self, z: ArrayLike):
        """Converts a depth expressed in meters in the corresponding
        index level.

        The returned value is an integer indicating the previous
        (*not* the nearest) depth in the z-levels.
        If `z` is above the first level, we return 0 anyway.

        Example:
            M = Mask(filename)
            k = M.get_depth_index(200.)
            M.zlevels[k]
        returns 192.60
        """
        return np.max(np.searchsorted(self._zlevels, z, side='right') - 1, 0)


class RegularMesh(Mesh, RegularGridDescriptor):
    def __init__(self, regular_grid: RegularGridDescriptor, zlevels: np.ndarray,
                 e3t: Optional[np.ndarray] = None):
        if not regular_grid.is_regular():
            raise ValueError(
                'regular_grid must be a regular grid.'
            )
        super().__init__(regular_grid, zlevels, e3t)

    @property
    def lon(self):
        # noinspection PyUnresolvedReferences
        return self._grid.lon

    @property
    def lat(self):
        # noinspection PyUnresolvedReferences
        return self._grid.lat
