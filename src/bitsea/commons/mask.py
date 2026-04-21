from os import PathLike
from pathlib import Path
from typing import Literal
from typing import Optional
from typing import Tuple
from typing import Union
from warnings import warn

import netCDF4
import numpy as np
import xarray as xr
from numpy.typing import ArrayLike
from scipy.sparse import dok_array

from bitsea.commons.bathymetry import Bathymetry
from bitsea.commons.geodistances import extend_from_average
from bitsea.commons.grid import Grid
from bitsea.commons.grid import MaskLayer
from bitsea.commons.grid import RegularGrid
from bitsea.commons.mesh import Mesh
from bitsea.components.component_mask_2d import ComponentMask2D
from bitsea.utilities.array_wrapper import BooleanArrayWrapper


# The fill value used for the missing values; this is a low precision number,
# so it is stable no matter what dtype is used for the numpy array
FILL_VALUE = np.float32(1e20)


class Mask(BooleanArrayWrapper, Mesh):
    """
    A `Mask` object represents the geometry of the domain in which the model
    operates.

    As a subclass of `Mesh`, it stores detailed information about cell geometry
    and the positions of cell centers. Additionally, it assigns a boolean value
    to each cell, indicating whether it contains data (water cell) or not (land
    cell).
    """

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

    def copy(self):
        """
        Returns a copy of the current object
        """
        return self.__class__(
            grid=self._grid.copy(),
            zlevels=self._zlevels.copy(),
            mask_array=self.as_array().copy(),
            allow_broadcast=False,
            e3t=self._e3t.copy(),
        )

    def get_sea_cells(self) -> np.ndarray:
        """
        Generates a 3D boolean array marking each water cell within the sea.

        For a standard `Mask` object, this output matches the mask content
        directly. However, for complex objects (e.g., those accounting for
        rivers), the result may differ, reflecting additional features.

        Returns:
            np.ndarray: A 3D array of boolean values, where `True` indicates
            a water cell located within the sea.
        """
        return self[:]

    def get_water_cells(self) -> np.ndarray:
        """
        Generates a 3D boolean array indicating the presence of water in each
        cell.

        For a standard `Mask` object, this output matches the full mask
        content. However, for specialized objects like a `SubMask`, which
        defines a subset of the original mask, only a portion of the water
        cells are marked as `True`. In this case, this method returns `True`
        also for the water cells that are outside the basin of the `SubMask`.

        Returns:
            np.ndarray: A 3D array of boolean values, where `True` indicates
            the presence of water in a cell.
        """
        return self[:]

    def convert_lon_lat_wetpoint_indices(
        self, *, lon: float, lat: float, max_radius: Optional[int] = 2
    ):
        """Converts longitude and latitude to the nearest water point index
        on the mask with maximum distance limit.

        Be aware that this function returns the indices following the "first
        lon then lat" convention. Therefore, they are in the opposite order
        respect to the 2d arrays returned by this object

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
        jp, ip = self.convert_lat_lon_to_indices(lon=lon, lat=lat)
        if self[0, jp, ip]:
            return ip, jp

        if max_radius is None:
            max_radius = max(self.shape[1], self.shape[2])

        # Candidate indices where we can look for a wet point
        # j_indices is the i indices of the mask that must be analyzed when
        # we look for the minimum (excluding the points that are too far).
        # Later we will do the same also for j. Moreover, ip_position is the
        # position of our first candidate inside the new array of indices that
        # we have just created. This is useful because, if we work on a small
        # portion of the domain, we need a way to go back to the original
        # indices. For example, for every 1D array
        # v[ip] = v[j_indices][ip_position]
        left_j_side = max(jp - max_radius, 0)
        right_j_side = min(jp + max_radius + 1, self.shape[1])
        j_indices = np.arange(left_j_side, right_j_side)
        jp_position = left_j_side

        left_i_side = max(ip - max_radius, 0)
        right_i_side = min(ip + max_radius + 1, self.shape[2])
        i_indices = np.arange(left_i_side, right_i_side)
        ip_position = left_i_side

        # We cut the mask around the point we have found
        i_slice = slice(i_indices[0], i_indices[-1] + 1)
        j_slice = slice(j_indices[0], j_indices[-1] + 1)
        cut_mask = self[0, j_slice, i_slice].copy()

        # We create a 2D array that saves the distances of these points
        # from (jp, ip) in grid units
        j_dist = j_indices - jp
        j_dist *= j_dist

        i_dist = i_indices - ip
        i_dist *= i_dist

        distances = np.sqrt(i_dist + j_dist[:, np.newaxis])

        # Remove all the points whose distance is greater than the max_radius
        cut_mask[distances > max_radius] = False

        # If there is no water in the slice, we return (jp, ip), but we
        # also raise a warning
        if not np.any(cut_mask):
            warn(
                f"Around point (lon={lon}, lat={lat}), whose closest grid "
                f"point is ({jp}, {ip}), there is no water point in a radius "
                f"of {max_radius} grid points"
            )
            return ip, jp

        # Outside water, we fix a distance that is bigger than any other
        # distance
        distances[~cut_mask] = max_radius * max_radius + 1

        # We get the index of the minimum value. We need to unravel because
        # argmin works on the flattened array
        local_min = np.unravel_index(
            np.argmin(distances, axis=None), distances.shape
        )

        # The previous indices are computed on the cut_mask. Here we go back
        # to the original indices
        return ip_position + int(local_min[1]), jp_position + int(local_min[0])

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
        Returns a 2d map that associates for each column the number of water
        cells that are present on that column.

        Returns:
            np.ndarray: A 2d numpy array of integers
        """
        return np.count_nonzero(self._data_array, axis=0)

    def rough_bathymetry(self) -> np.ndarray:
        """
        Calculates the bathymetry used by the model
        It does not take into account e3t

        Returns:
            np.ndarray: a 2d numpy array of floats
        """
        n_cells_per_column = self.bathymetry_in_cells()
        zlevels = np.concatenate((np.array([0]), self.zlevels))
        return zlevels[n_cells_per_column]

    def bathymetry(self) -> np.ndarray:
        """
        Calculates the bathymetry used by the model
        Best evaluation, since it takes into account e3t.

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
        level_mask = self.mask_at_level(depth)

        # This is the value that we will return if there are no water cell
        # at the requested depth
        zero_array = np.zeros((0,), dtype=self.xlevels.dtype)

        # Ensure that there is at least one water cell
        if np.all(~level_mask):
            warn(
                f"At depth {depth}, there are no water points in the "
                "current mask"
            )
            return zero_array, zero_array.copy()

        # Find all the connected components of the sea
        component_mask = ComponentMask2D(level_mask)
        n_components = component_mask.n_components

        # Now we concatenate all the boundaries of each component. The output
        # will be a list of indices; we will later transform them into WSG84
        # coordinates
        boundary_vertex_indices = []
        for component in range(n_components):
            boundary_vertex_indices.extend(
                component_mask.get_component_boundary_curves(component)
            )

        # We must compute the position of the vertices of the cells; we do that
        # by imposing that the coordinates of the center of the cells are the
        # middle points of the coordinates of the vertices
        vertex_lat_coordinates = extend_from_average(
            extend_from_average(self.ylevels, axis=0), axis=1
        )
        vertex_lon_coordinates = extend_from_average(
            extend_from_average(self.xlevels, axis=0), axis=1
        )

        # Now we transform the paths that we have computed previously into
        # paths of coordinates. We also discard the paths that are shorter
        # than `min_line_length` and we save the length of the longest path
        longest_path = 0
        coord_paths = []
        for path in boundary_vertex_indices:
            longest_path = max(longest_path, len(path))
            if len(path) < min_line_length:
                continue
            lat_indices = tuple(p[0] for p in path)
            lon_indices = tuple(p[1] for p in path)
            path_lats = vertex_lat_coordinates[lat_indices, lon_indices]
            path_lons = vertex_lon_coordinates[lat_indices, lon_indices]

            coord_paths.append((path_lats, path_lons))

        if len(coord_paths) == 0:
            # The filter on the min_line_length has removed all the points
            # We return two empty arrays
            warn(
                f"min_line_length is set to {min_line_length}, but at depth "
                f"{depth}m the longest boundary path is {longest_path}; "
                f"the output of the method `coastline` will be empty"
            )
            return zero_array, zero_array.copy()

        # Now we have the paths; we want to concatenate them inserting a nan
        # between two different paths. In this way, it will be easy to plot
        # the results (because the nan will split the line of matplotlib).
        lat_values = (p[0] for p in coord_paths)
        lon_values = (p[1] for p in coord_paths)

        # We prepare the array that will be used to split: it contains just
        # one nan
        nan_split = np.full((1,), fill_value=np.nan)

        lat_separated_values = []
        for lat_path in lat_values:
            lat_separated_values.append(lat_path)
            lat_separated_values.append(nan_split)
        lat_separated_values = tuple(lat_separated_values[:-1])

        lon_separated_values = []
        for lon_path in lon_values:
            lon_separated_values.append(lon_path)
            lon_separated_values.append(nan_split)
        lon_separated_values = tuple(lon_separated_values[:-1])

        # Finally, we can concatenate our results
        lat_values = np.concatenate(lat_separated_values)
        lon_values = np.concatenate(lon_separated_values)

        return lon_values, lat_values

    def cut_at_level(self, index: int):
        """
        Retrieves a `MaskLayer` at the specified depth level.

        Arguments:
            index (int): The depth level index.

        Returns:
            MaskLayer: The `MaskLayer` corresponding to the specified depth
              index.
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

    def column_side_area(
        self,
        ji: int,
        jj: int,
        side: Literal["N", "S", "W", "E"],
        n_vertical_cells: Optional[int] = None,
    ) -> float:
        """
        Calculates the lateral area of a water column with the specified depth.

        Args:
            ji (int): Longitudinal index.
            jj (int): Latitudinal index.
            side (Literal["N", "S", "W", "E"]): The side of the column for
              which to compute the area, specified by the first letter of a
              cardinal direction.
            n_vertical_cells (Optional[int]): The number of cells in the
              column. If it is `None`, the column goes from the surface to
              the bottom layer

        Returns:
            float: The lateral area of the specified column side.
        """
        if n_vertical_cells is None:
            n_vertical_cells = np.count_nonzero(self[:, jj, ji], axis=0)

        return super().column_side_area(ji, jj, side, n_vertical_cells)

    def to_xarray(self) -> xr.Dataset:
        """
        Converts the `Mask` object to an `xarray.Dataset`. The new dataset
        follows the name conventions of CMEMS for the coordinates: the
        latitude is called "latitude" (instead of "lat"), the "longitude" is
        called "longitude", and the depth is called "depth".

        Returns:
            An Xarray Dataset
        """
        if self.grid.is_regular():
            return xr.Dataset(
                {
                    "tmask": (("depth", "latitude", "longitude"), self),
                },
                coords={
                    "latitude": (("latitude",), self.lat),
                    "longitude": (("longitude",), self.lon),
                    "depth": (("depth",), self.zlevels),
                },
            )
        return xr.Dataset(
            {
                "tmask": (("depth", "y", "x"), self),
            },
            coords={
                "latitude": (("y", "x"), self.ylevels),
                "longitude": (("y", "x"), self.xlevels),
                "depth": (("depth",), self.zlevels),
            },
        )

    @staticmethod
    def from_xarray(dataset: xr.Dataset) -> "Mask":
        """
        Create a Mask object from an Xarray Dataset.

        This function is the inverse of `Mask.to_xarray`. It takes an Xarray
        generated from a *regular* Mask (i.e., one that has a regular grid)
        using its method `to_xarray` and returns a Mask object.

        This function does not support Datasets generated from a `Mask` that
        is not regular.

        Args:
            dataset (xr.Dataset): The input xarray Dataset containing latitude,
                                  longitude, depth, and mask data.

        Returns:
            Mask: A Mask object constructed using the provided xarray Dataset.
        """
        latitude = dataset.latitude.values
        longitude = dataset.longitude.values
        grid = RegularGrid(lat=latitude, lon=longitude)

        depth = dataset.depth.values
        return Mask(grid=grid, zlevels=depth, mask_array=dataset.tmask.values)

    @classmethod
    def from_file_pointer(
        cls,
        file_pointer: netCDF4.Dataset,
        *,
        zlevels_var_name: str = "nav_lev",
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
        e3t_var_name: Optional[str] = None,
        mask_var_name: str = "tmask",
        read_e3t: bool = True,
    ):
        """
        Creates a new `Mask` object from a netCDF file.

        Args:
            file_pointer (netCDF4.Dataset): The netCDF file containing mask
              data.
            zlevels_var_name (str): The name of the variable representing
              depth levels (zlevels) in the file.
            ylevels_var_name (str): The name of the variable representing
              latitude coordinates (ylevels) in the file.
            xlevels_var_name (str): The name of the variable representing
              longitude coordinates (xlevels) in the file.
            e3t_var_name (Optional[str]): The name of the variable representing
              the vertical size of the levels in the file.
            mask_var_name (str): The name of the variable representing the
              mask in the file.
            read_e3t (bool): If `True`, reads the `e3t` variable from the
              file; if `False`, computes `e3t` from `zlevels`.

        Returns:
            Mask: A new `Mask` object created from the specified netCDF data.
        """
        mesh = Mesh.from_file_pointer(
            file_pointer,
            zlevels_var_name=zlevels_var_name,
            ylevels_var_name=ylevels_var_name,
            xlevels_var_name=xlevels_var_name,
            e3t_var_name=e3t_var_name,
            read_e3t=read_e3t,
        )

        mask_array = np.asarray(
            file_pointer.variables[mask_var_name][:], dtype=bool
        )

        if mask_array.ndim == 4:
            mask_array = mask_array[0, :, :, :]

        if read_e3t:
            e3t = mesh.e3t
        else:
            e3t = None

        return cls(mesh.grid, mesh.zlevels, mask_array, e3t)

    @classmethod
    def from_file(
        cls,
        file_path: PathLike,
        *,
        zlevels_var_name: str = "nav_lev",
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
        e3t_var_name: Optional[str] = None,
        mask_var_name: str = "tmask",
        read_e3t: bool = True,
    ):
        """
        Creates a new `Mask` object from a netCDF file.

        This method is similar to `from_file_pointer`, but instead of
        accepting a netCDF `Dataset` object, it takes the file path as input.

        Args:
            file_path (PathLike): The path to the netCDF file containing
              the mask.
            zlevels_var_name (str): The name of the variable representing
              depth levels (zlevels) in the file.
            ylevels_var_name (str): The name of the variable representing
              latitude coordinates (ylevels) in the file.
            xlevels_var_name (str): The name of the variable representing
              longitude coordinates (xlevels) in the file.
            e3t_var_name (Optional[str]): The name of the variable representing
              the vertical size of the levels in the file.
            mask_var_name (str): The name of the variable representing the
              mask in the file.
            read_e3t (bool): If `True`, reads the `e3t` variable from the
              file; if `False`, computes `e3t` from `zlevels`.

        Returns:
            Mask: A new `Mask` object created from the specified netCDF file.
        """
        with netCDF4.Dataset(file_path, "r") as f:
            return cls.from_file_pointer(
                f,
                zlevels_var_name=zlevels_var_name,
                ylevels_var_name=ylevels_var_name,
                xlevels_var_name=xlevels_var_name,
                e3t_var_name=e3t_var_name,
                mask_var_name=mask_var_name,
                read_e3t=read_e3t,
            )

    @classmethod
    def from_mer_file(
        cls,
        file_path: PathLike,
    ):
        with xr.open_dataset(file_path) as ds:
            latitude = ds.latitude.values
            longitude = ds.longitude.values
            zlevels = ds.depth.values
            mask_array = ds.tmask.values == 1

            if "e3t" in ds.variables:
                e3t = ds.e3t.values
            else:
                e3t = None

        grid = RegularGrid(lat=latitude, lon=longitude)
        mask = cls(grid=grid, zlevels=zlevels, mask_array=mask_array, e3t=e3t)
        return mask

    def _add_attributes_on_file(self, file_pointer):
        """
        Called by the `save_as_netcdf` method.

        This method does nothing by default, but subclasses of `Mask` can
        override it to add custom attributes to the netCDF file generated by
        `save_as_netcdf`.

        Args:
            file_pointer (netCDF4.Dataset): The netCDF dataset where attributes
             will be saved.
        """
        pass

    def _save_as_netcdf(self, file_pointer: netCDF4.Dataset, compression_level):
        """
        Internal method to save the current `Mask` object to a netCDF file.

        This is the internal implementation of the `save_as_netcdf` method.
        Unlike the public method, which opens the file internally, this method
        operates on an already-open file pointer.

        Subclasses inheriting from `Mask` can override this method to include
        additional information in the output file.

        Args:
            file_pointer (netCDF4.Dataset): An open netCDF Dataset object where
                the mask data will be written.
            compression_level (int): Compression level for the saved variables.
                Use 0 for no compression, or an integer between 1 and 9, where
                1 provides the fastest compression and 9 provides the highest
                compression ratio. The default value is 3, offering a good
                balance between compression speed and file size.
        """
        var_kwargs = {}
        if compression_level > 0:
            var_kwargs["zlib"] = True
            var_kwargs["complevel"] = compression_level

        # Add the spatial dimensions
        file_pointer.createDimension("t", 1)
        file_pointer.createDimension("z", self.shape[0])
        file_pointer.createDimension("y", self.shape[1])
        file_pointer.createDimension("x", self.shape[2])

        # Add nav_lev data
        nav_lev = file_pointer.createVariable(
            "nav_lev", "f4", ("z",), **var_kwargs
        )
        nav_lev[:] = self.zlevels

        # Create the nav_lat NetCDF variable
        nav_lat = file_pointer.createVariable(
            "nav_lat", "f4", ("y", "x"), **var_kwargs
        )
        nav_lat[:, :] = self.ylevels

        # Create the nav_lon NetCDF variable
        nav_lon = file_pointer.createVariable(
            "nav_lon", "f4", ("y", "x"), **var_kwargs
        )
        nav_lon[:, :] = self.xlevels

        e1t = file_pointer.createVariable("e1t", "f4", ("y", "x"), **var_kwargs)
        e1t[:] = self.e1t

        e2t = file_pointer.createVariable("e2t", "f4", ("y", "x"), **var_kwargs)
        e2t[:] = self.e2t

        e3t = file_pointer.createVariable(
            "e3t", "f4", ("t", "z", "y", "x"), **var_kwargs
        )
        e3t[:] = self.e3t

        # Create a variable to hold the data
        mask = file_pointer.createVariable(
            "tmask", "u1", ("z", "y", "x"), **var_kwargs
        )

        # Assign the data to the NetCDF variable
        mask[:, :, :] = self[:, :, :]

        self._add_attributes_on_file(file_pointer)

    def save_as_netcdf(
        self,
        file_path: Union[PathLike, str],
        compression_level: Literal[0, 1, 2, 3, 4, 5, 6, 7, 8, 9] = 3,
    ):
        """
        Saves the Mask into a netCDF file.

        Saves the current `Mask` object to a netCDF file, which can later
        be loaded using the `from_file` method.

        Args:
            file_path: The path to the netCDF file where the `Mask` object
                will be saved.
            compression_level: The compression level to use when saving the
                variables. Valid values are 0 (no compression), and integers
                between 1 and 9, where 1 is the fastest compression and 9 is
                the slowest. The default value is 3, which is a reasonable
                compromise between speed and compression.
        """
        file_path = Path(file_path)
        if compression_level < 0 or compression_level > 9:
            raise ValueError(
                "compression_level must be an integer between 0 and 9; "
                f"received {compression_level}."
            )

        with netCDF4.Dataset(file_path, "w") as netCDF_out:
            self._save_as_netcdf(
                netCDF_out, compression_level=compression_level
            )


class MaskBathymetry(Bathymetry):
    """
    This class is a bathymetry generated starting from a mask, i.e., it
    returns the z-coordinate of the bottom face of the deepest cell of the
    column that contains the point (lon, lat).
    """

    def __init__(self, mask: Mask):
        self._mask = mask
        self._bathymetry_data = mask.bathymetry()

        # Fix the bathymetry of the land cells to 0 (to be coherent with the
        # behavior of the bathymetry classes). Otherwise, if we let the land
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
            indices = self._mask.convert_lat_lon_to_indices(
                lon=lon[inside_domain], lat=lat[inside_domain]
            )
            output[inside_domain] = self._bathymetry_data[indices]

        if not return_float:
            return output

        return float(output.squeeze().item())


class MaskWithRivers(Mask):
    """
    A subclass of Mask that distinguishes river cells from sea and land cells.

    This class extends the base Mask functionality by allowing certain cells to
    be marked as rivers, in addition to the standard sea and land
    classification.
    River cells are handled separately from sea cells, which is useful for
    models or analyses that require to avoid river cells.
    """

    def __init__(
        self,
        grid: Grid,
        zlevels: ArrayLike,
        mask_array: ArrayLike,
        river_positions: ArrayLike,
        allow_broadcast: bool = False,
        e3t: Optional[np.ndarray] = None,
    ):
        """Initialize a MaskWithRivers object.

        Args:
            grid: The grid object defining the spatial domain.
            zlevels: The vertical levels of the mask (e.g., depth layers).
            mask_array: A boolean array indicating the presence of
              water (True) or land (False) at each grid cell and depth level.
              Shape should be (len(zlevels), grid.shape[0], grid.shape[1]).
              The mask shall be `True` also on the cells occupied by rivers.
            river_positions: An array or sparse matrix indicating the
              positions of river cells in the grid. Typically, this should
              be a 2D array (grid.shape[0], grid.shape[1]) where nonzero
              entries indicate river cells and their values may represent
              river indices or IDs.
            allow_broadcast (bool): If True, allows broadcasting of mask_array
              to match the required shape. Default is False.
            e3t (Optional[np.ndarray]): Optional array specifying the
              thickness of each vertical layer.

        Note:
            River cells are removed from the mask and stored separately for
            further processing.
        """
        rivers = np.asarray(river_positions)

        mask_dim = len(zlevels), grid.shape[0], grid.shape[1]

        try:
            mask_data = np.asarray(mask_array, dtype=bool, copy=False)
            input_copied = False
        except TypeError:
            # Old versions of Numpy do not support the "copy=False" syntax;
            # In this case, we call the function without specifying the value
            # of the "copy" argument, and we must assume that the data has not
            # been copied (for safety)
            mask_data = np.asarray(mask_array, dtype=bool)
            input_copied = False
        except ValueError:
            # In this case, it has been impossible to avoid the copy. We run
            # again the same command, but this time we force the copy
            mask_data = np.asarray(mask_array, dtype=bool)
            input_copied = True

        # If mask_data is too small to cover all the mask, we try to broadcast
        # it. This enforces a copy
        if allow_broadcast and mask_data.shape != mask_dim:
            mask_data = np.copy(np.broadcast_to(mask_array, mask_dim))
            input_copied = True

        # Ensure that mask_array is writeable
        if not input_copied:
            mask_data = np.copy(mask_data)

        # Here we save the cells that belong to each river into a list of
        # sparse arrays. We must use a list of 2D array because there are no
        # 3D sparse arrays implemented inside scipy. So the index of the list
        # is the depth
        self._river_cells = []
        rivers = dok_array(rivers)
        for depth_index in range(mask_data.shape[0]):
            found_cells = False
            current_river_cells = dok_array(mask_data.shape[1:], dtype=int)
            for (i, j), value in rivers.items():
                if not mask_data[depth_index, i, j]:
                    continue
                current_river_cells[i, j] = value
                found_cells = True

            if not found_cells:
                break
            self._river_cells.append(current_river_cells)

        # Remove rivers from mask_array
        for i, j in zip(*rivers.nonzero()):
            mask_data[:, i, j] = False

        super().__init__(
            grid=grid,
            zlevels=zlevels,
            mask_array=mask_data,
            allow_broadcast=False,
            e3t=e3t,
        )

    def copy(self):
        """
        Return a copy of this MaskWithRivers, preserving river information.
        """

        return self.__class__(
            grid=self.grid,
            zlevels=self.zlevels,
            river_positions=self._river_cells[0].todense(),
            mask_array=self.get_water_cells(),
            allow_broadcast=False,
            e3t=getattr(self, "e3t", None),
        )

    def get_water_cells(self) -> np.ndarray:
        """Return a mask of all water cells: both sea and river cells.

        Unlike the parent class's mask (which excludes river cells), this
        method returns the combined mask of sea cells and river cells. This is
        useful when you need to treat rivers as part of the water domain.

        Returns:
            A boolean array of the same shape as the mask, where
            True indicates a water cell (either sea or river) and False
            indicates a land cell.
        """
        output_data = np.copy(self[:])
        for depth_index, river_cells in enumerate(self._river_cells):
            lat_indices, lon_indices = river_cells.nonzero()
            output_data[depth_index, lat_indices, lon_indices] = True
        return output_data

    def get_river_ids(self) -> set:
        """Return a set of unique river identifiers present in the mask.

        This method extracts the unique river IDs from the river positions
        stored in the mask. It assumes that river cells are labeled with
        integer values representing different rivers.

        Returns:
            set: A set of unique river identifiers found in the mask.
        """
        if len(self._river_cells) == 0:
            return set()

        # Get the unique river IDs from the first depth layer
        unique_rivers = set(self._river_cells[0].values())

        return unique_rivers

    def n_rivers(self) -> int:
        """Return the number of distinct rivers represented in the mask.

        This method counts the number of unique river identifiers present in
        the river positions. It assumes that river cells are labeled with
        integer values representing different rivers.

        Returns:
            int: The number of distinct rivers in the mask.
        """
        return len(self.get_river_ids())

    @classmethod
    def from_file_pointer(
        cls,
        file_pointer: netCDF4.Dataset,
        *,
        zlevels_var_name: str = "nav_lev",
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
        e3t_var_name: Optional[str] = None,
        mask_var_name: str = "tmask",
        rivers_var_name: str = "rivers",
        read_e3t: bool = True,
    ):
        raw_mask = Mask.from_file_pointer(
            file_pointer,
            zlevels_var_name=zlevels_var_name,
            ylevels_var_name=ylevels_var_name,
            xlevels_var_name=xlevels_var_name,
            e3t_var_name=e3t_var_name,
            mask_var_name=mask_var_name,
            read_e3t=read_e3t,
        )
        rivers = np.asarray(
            file_pointer.variables[rivers_var_name][:], dtype=int
        )
        return cls(
            grid=raw_mask.grid,
            zlevels=raw_mask.zlevels,
            mask_array=raw_mask.as_mutable_array(),
            river_positions=rivers,
            allow_broadcast=False,
            e3t=raw_mask.e3t,
        )

    @classmethod
    def from_file(
        cls,
        file_path: PathLike,
        *,
        zlevels_var_name: str = "nav_lev",
        ylevels_var_name: str = "nav_lat",
        xlevels_var_name: str = "nav_lon",
        e3t_var_name: Optional[str] = None,
        mask_var_name: str = "tmask",
        rivers_var_name: str = "rivers",
        read_e3t: bool = True,
    ):
        with netCDF4.Dataset(file_path, "r") as f:
            return cls.from_file_pointer(
                f,
                zlevels_var_name=zlevels_var_name,
                ylevels_var_name=ylevels_var_name,
                xlevels_var_name=xlevels_var_name,
                e3t_var_name=e3t_var_name,
                mask_var_name=mask_var_name,
                rivers_var_name=rivers_var_name,
                read_e3t=read_e3t,
            )

    def as_standard_mask(self, include_rivers: bool = False) -> Mask:
        """
        Converts the MaskWithRivers object to a standard Mask object.

        This method creates a new standard Mask object from the current
        MaskWithRivers instance. By default, river cells are excluded from the
        resulting mask (treated as land). If `include_rivers` is set to True,
        river cells are included as water cells in the output mask.

        Args:
            include_rivers (bool): If True, includes river cells as water cells
                in the output mask. If False (default), river cells are treated
                as land cells.

        Returns:
            Mask: A standard Mask object with the same grid and vertical levels
                as the current MaskWithRivers instance.
        """
        if include_rivers:
            data_array = self.get_water_cells()
        else:
            data_array = self.as_array()

        return Mask(
            grid=self.grid,
            zlevels=self.zlevels,
            mask_array=data_array,
            e3t=self.e3t,
        )

    def save_as_netcdf(
        self,
        file_path: Union[PathLike, str],
        compression_level: Literal[0, 1, 2, 3, 4, 5, 6, 7, 8, 9] = 3,
    ):
        file_path = Path(file_path)
        if compression_level < 0 or compression_level > 9:
            raise ValueError(
                "compression_level must be an integer between 0 and 9; "
                f"received {compression_level}."
            )

        var_kwargs = {}
        if compression_level > 0:
            var_kwargs["zlib"] = True
            var_kwargs["complevel"] = compression_level

        with netCDF4.Dataset(file_path, "w") as netCDF_out:
            self.as_standard_mask(include_rivers=True)._save_as_netcdf(
                netCDF_out, compression_level=compression_level
            )
            var_kwargs["fill_value"] = 0
            rivers_var = netCDF_out.createVariable(
                "rivers", np.int32, ("y", "x"), **var_kwargs
            )
            if self.n_rivers() > 0:
                rivers_var[:] = self._river_cells[0].toarray()
            else:
                rivers_var[:] = 0
