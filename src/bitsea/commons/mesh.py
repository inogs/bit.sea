from typing import Optional, Tuple, Union

import numpy as np

from bitsea.commons.grid import GridDescriptor
from bitsea.commons.grid import extend_from_average


class Mesh(GridDescriptor):
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

        if e3t is not None:
            expected_shape = (self._zlevels.shape[0],) + self._grid.shape
            self._e3t = np.asarray(e3t).view()
            if self._e3t.shape != expected_shape:
                raise ValueError(
                    f'e3t is expected to be {expected_shape}, but got '
                    f'an array with shape {self._e3t.shape}'
                )
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
    def xlevels(self) -> np.ndarray:
        return self._grid.xlevels

    @property
    def ylevels(self) -> np.ndarray:
        return self._grid.ylevels

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

    def e3t(self) -> np.ndarray:
        return self._e3t

    def convert_lon_lat_to_indices(self, *, lon: Union[float, np.ndarray],
                                   lat: Union[float, np.ndarray]) -> Tuple:
        return self._grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    def convert_i_j_to_lon_lat(self, i: Union[int, np.ndarray],
                               j: Union[int, np.ndarray]) -> Tuple:
        return self._grid.convert_i_j_to_lon_lat(i=i, j=j)
