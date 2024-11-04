from os import PathLike
from typing import Optional
from typing import Tuple
from warnings import warn

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from numpy.typing import ArrayLike

from bitsea.commons.bathymetry import Bathymetry
from bitsea.commons.grid import Grid
from bitsea.commons.grid import MaskLayer
from bitsea.commons.mesh import Mesh
from bitsea.commons.mesh import RegularMesh
from bitsea.utilities.array_wrapper import BooleanArrayWrapper


# The fill value used for the missing values
FILL_VALUE = 1e20


class Mask(BooleanArrayWrapper, Mesh):
    def __init__(
        self,
        grid: Grid,
        zlevels: ArrayLike,
        mask_array: ArrayLike,
        allow_broadcast: bool = False,
        e3t: Optional[np.ndarray] = None,
    ):
        mask_array = np.asarray(mask_array)
        zlevels = np.asarray(zlevels)

        if zlevels.ndim != 1:
            raise ValueError(
                f"zlevels must be a 1D array; its current shape is "
                f"{zlevels.shape}"
            )

        expected_shape = (zlevels.shape[0],) + grid.shape
        if mask_array.shape != expected_shape:
            if allow_broadcast:
                mask_array = np.broadcast_to(mask_array, expected_shape)
            else:
                raise ValueError(
                    f"mask_array has shape {mask_array.shape} while the "
                    f"expected shape was {expected_shape}"
                )

        super().__init__(
            grid=grid, zlevels=zlevels, e3t=e3t, wrapped_data=mask_array
        )

    @property
    def mask(self):
        warn(
            "This method is deprecated. Use the object itself instead",
            DeprecationWarning,
        )
        return self

    def get_sea_cells(self):
        return self[:]

    def get_water_cells(self):
        return self[:]

    def convert_lon_lat_wetpoint_indices(
        self, *, lon: float, lat: float, max_radius: Optional[int] = 2
    ):
        """Converts longitude and latitude to the nearest water point index
        on the mask with maximum distance limit

        Args:
            lon (float): Longitude in degrees.
            lat (float): Latitude in degrees.
            max_radius (Optional[int]): Maximum distance where the water point
              is searched (in grid units, integer, default: 2)

        Returns:
            a tuple of numbers, the first one is the longitude index and
            the other one is the latitude index.
        """

        # Indexes of the input lon, lat
        jp, ip = self.convert_lon_lat_to_indices(lon=lon, lat=lat)
        if self[0, jp, ip]:
            return jp, ip

        if max_radius is None:
            max_radius = max(self.shape[1], self.shape[2])

        # Candidate indices where we can look for a wet point
        left_i_side = max(ip - max_radius, 0)
        right_i_side = min(ip + max_radius + 1, self.shape[2])
        i_indices = np.arange(left_i_side, right_i_side)
        ip_position = left_i_side

        left_j_side = max(jp - max_radius, 0)
        right_j_side = min(jp + max_radius + 1, self.shape[1])
        j_indices = np.arange(left_j_side, right_j_side)
        jp_position = left_j_side

        # We cut the mask around the point we have found
        j_slice = slice(j_indices[0], j_indices[-1] + 1)
        i_slice = slice(i_indices[0], i_indices[-1] + 1)
        cut_mask = self[0, j_slice, i_slice].copy()

        # We create a 2D array that saves the distances of these points
        # from (jp, ip) in grid units
        i_dist = i_indices - ip
        i_dist *= i_dist

        j_dist = j_indices - jp
        j_dist *= j_dist

        distances = np.sqrt(i_dist + j_dist[:, np.newaxis])

        # Remove all the points whose distance is greater than the max_radius
        cut_mask[distances > max_radius] = False

        # If there is no water in the slice, we return (jp, ip) but we
        # also raise a warning
        if not np.any(cut_mask):
            warn(
                f"Around point (lon={lon}, lat={lat}), whose closest grid "
                f"point is ({jp}, {ip}), there is no water point in a radius "
                f"of {max_radius} grid points"
            )
            return jp, ip

        # Outside water, we fix a distance that is bigger than any other
        # distance
        distances[~cut_mask] = max_radius * max_radius + 1

        local_min = np.unravel_index(
            np.argmin(distances, axis=None), distances.shape
        )

        return jp_position + int(local_min[0]), ip_position + int(local_min[1])

    def mask_at_level(self, z: float) -> np.ndarray:
        """
        Returns a 2d slice of the tmask for the depth (m) provided as argument.

        When z is not a point in the discretization, jk_m is selected as the
        point immediately before. This depth level should not be included in
        the mask:

        (jk_m-1)   FFFFFFFFFFFFFFFF
        (jk_m  )   FFFFFFFFFFFFFFFF
        z----------------------
        (jk_m +1)  TTTTTTTTTTTTTTTT
        (jk_m +2)  TTTTTTTTTTTTTTTT
        """
        if z < self.zlevels[0]:
            return self[0, :, :]

        jk_m = self.get_depth_index(z)

        # if we are beyond the last level, return False everywhere
        if jk_m == self.zlevels.shape[0] - 1:
            return np.zeros_like(self)

        return self[jk_m + 1, :, :]

    def bathymetry_in_cells(self) -> np.ndarray:
        """
        Returns a 2d map that for each columns associates the number of water
        cells that are present on that column.

        Returns:
            A 2d numpy array of integers
        """
        return np.count_nonzero(self._data_array, axis=0)

    def rough_bathymetry(self):
        """
        Calculates the bathymetry used by the model
        It does not takes in account e3t

        Returns:
        * bathy * a 2d numpy array of floats
        """
        n_cells_per_column = self.bathymetry_in_cells()
        zlevels = np.concatenate((np.array([0]), self.zlevels))
        return zlevels[n_cells_per_column]

    def bathymetry(self):
        """
        Calculates the bathymetry used by the model
        Best evaluation, since it takes in account e3t.

        Returns:
          a 2d numpy array of floats
        """
        bathymetry = np.sum(self.e3t, axis=0, where=self[:])
        bathymetry[~self[0, :, :]] = FILL_VALUE
        return bathymetry

    def coastline(
        self, depth: float, min_line_length: int = 30
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculates a mesh-dependent coastline at a chosen depth.

        Arguments :
          depth (float): depths expressed in meters
          min_line_length (int): integer indicating the minimum number of
            points of each coastline, it is a threshold to avoid a lot of
            small islands.

        Returns:
          (x,y): A tuple with two numpy 1d arrays, containing NaNs to separate
            the lines, in order to be easily plotted.
        """
        level_mask = self.mask_at_level(depth).astype(np.float64)

        mask_contour = plt.contour(
            self.xlevels, self.ylevels, level_mask, levels=[float(0.5)]
        )

        path_list = mask_contour.collections[0].get_paths()

        point_coords = np.zeros((0, 2))
        nan = np.full((1, 2), np.nan, dtype=np.float32)

        for p in path_list:
            v = p.vertices
            n_points, _ = v.shape
            if n_points > min_line_length:
                point_coords = np.concatenate((point_coords, v, nan), axis=0)

        return point_coords[:-1, 0], point_coords[:-1, 1]

    def cut_at_level(self, index):
        """
        Arguments:
        * index * integer, depth index

        Returns copy of the mask object, reduced on a single horizontal slice,
        the one of the provided depth index.
        """
        return MaskLayer.from_grid(
            self.grid,
            depth=self.zlevels[index],
            thickness=self.dz[index],
            mask=self[index, :],
        )

    def CellArea(self, ji: int, jj: int, side, n_vertical_cells: int):
        warn(
            'This method has been renamed "column_side_are".',
            DeprecationWarning,
        )
        return self.column_side_area(ji, jj, side, n_vertical_cells)

    def column_side_area(self, ji: int, jj: int, side, n_vertical_cells: int):
        """
        Returns the lateral area of watercolumn of n Vertical cells.
        """
        # ji,jj = self.convert_lon_lat_to_indices(lon, lat)
        if side in ["E", "W"]:
            return self.e2t[jj, ji] * self.e3t[:n_vertical_cells, jj, ji].sum()
        elif side in ["N", "S"]:
            return self.e1t[jj, ji] * self.e3t[:n_vertical_cells, jj, ji].sum()
        else:
            raise ValueError(f'Invalid side: "{side}"')

    @classmethod
    def from_file_pointer(
        cls,
        file_pointer: netCDF4.Dataset,
        zlevels_var_name: str = "nav_lev",
        mask_var_name: str = "tmask",
        read_e3t: bool = True,
    ):
        mesh = Mesh.from_file_pointer(file_pointer, zlevels_var_name, read_e3t)

        mask_array = np.asarray(
            file_pointer.variables[mask_var_name][:], dtype=bool
        )

        if mask_array.ndim == 4:
            mask_array = mask_array[0, :, :, :]

        if read_e3t:
            e3t = mesh.e3t
        else:
            e3t = None

        if mesh.is_regular():
            return RegularMask(mesh.grid, mesh.zlevels, mask_array, e3t)

        return cls(mesh.grid, mesh.zlevels, mask_array, e3t)

    @classmethod
    def from_file(
        cls,
        file_path: PathLike,
        zlevels_var_name: str = "nav_lev",
        mask_var_name: str = "tmask",
        read_e3t: bool = True,
    ):
        with netCDF4.Dataset(file_path, "r") as f:
            return cls.from_file_pointer(
                f, zlevels_var_name, mask_var_name, read_e3t
            )


class RegularMask(Mask, RegularMesh):
    def __init__(
        self,
        grid: Grid,
        zlevels: ArrayLike,
        mask_array: ArrayLike,
        e3t: Optional[np.ndarray] = None,
    ):
        if not grid.is_regular():
            raise ValueError("Grid must be regular")

        super().__init__(
            grid=grid, zlevels=zlevels, mask_array=mask_array, e3t=e3t
        )


class MaskBathymetry(Bathymetry):
    """
    This class is a bathymetry,  generated starting from a mask, i.e., it
    returns the z-coordinate of the bottom face of the deepest cell of the
    column that contains the point (lon, lat).
    """

    def __init__(self, mask: Mask):
        self._mask = mask
        self._bathymetry_data = mask.bathymetry()

        # Fix the bathymetry of the land cells to 0 (to be coherent with the
        # behaviour of the bathymetry classes). Otherwise, if we let the land
        # points to have bathymetry = 1e20, they will be in every
        # BathymetricBasin
        self._bathymetry_data[np.logical_not(self._mask[0, :, :])] = 0

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, repr(self._mask))

    def is_inside_domain(self, lon, lat):
        return self._mask.is_inside_domain(
            lon=lon, lat=lat, raise_if_outside=False
        )

    def __call__(self, lon, lat):
        return_float = True
        if hasattr(lon, "__len__") or hasattr(lat, "__len__"):
            return_float = False

        lon, lat = np.broadcast_arrays(lon, lat)
        inside_domain = self.is_inside_domain(lon, lat)

        output = np.zeros_like(lon, dtype=np.float32)
        if np.any(inside_domain):
            indices = self._mask.convert_lon_lat_to_indices(
                lon=lon[inside_domain], lat=lat[inside_domain]
            )
            output[inside_domain] = self._bathymetry_data[indices]

        if not return_float:
            return output

        return float(output.squeeze().item())
