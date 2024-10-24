from abc import ABC
from typing import Optional
from typing import Tuple
from warnings import warn

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from scipy import spatial

from bitsea.commons.bathymetry import Bathymetry
from bitsea.commons.grid import Grid
from bitsea.commons.mesh import Mesh
from bitsea.commons.mesh import RegularMesh
from bitsea.commons.utils import search_closest_sorted
from bitsea.utilities.array_wrapper import BooleanArrayWrapper


class Mask(Mesh, BooleanArrayWrapper):
    def __init__(
        self,
        grid: Grid,
        zlevels: np.ndarray,
        mask_array: np.ndarray,
        e3t: Optional[np.ndarray] = None,
    ):
        super().__init__(grid, zlevels, e3t)

        if mask_array.shape != self.shape:
            mask_array = np.broadcast_to(mask_array, self.shape)

        self._mask_array = mask_array.view(dtype=bool)

    @property
    def mask(self):
        warn(
            "This method is deprecated. Use the object itself instead",
            DeprecationWarning,
        )
        return self

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
            return ip, jp

        if max_radius is None:
            max_radius = max(self.shape[1], self.shape[2])

        # Candidate indices where we can look for a wetpoint
        i_indices = np.arange(
            max(ip - max_radius, 0), min(ip + max_radius + 1, self.shape[2])
        )

        j_indices = np.arange(
            max(jp - max_radius, 0), min(jp + max_radius + 1, self.shape[2])
        )

        # We cut the mask around the point we have found
        j_slice = slice(j_indices[0], j_indices[-1] + 1)
        i_slice = slice(i_indices[0], i_indices[-1] + 1)
        cut_mask = self[0, j_slice, i_slice]

        # If there is no water in this position, we return (jp, ip) but we
        # also raise a warning
        if not np.any(cut_mask):
            warn(
                f"Around point (lon={lon}, lat={lat}), whose closest grid "
                f"point is ({jp}, {ip}), there is no water point in a radius "
                f"of {max_radius} grid points"
            )
            return jp, ip

        # We create a 2D array that saves the distances of these points
        # from (jp, ip) in grid units
        i_dist = i_indices - ip
        i_dist *= i_dist

        j_dist = j_indices - jp
        j_dist *= j_dist

        distances = np.sqrt(i_dist + j_dist[:, np.newaxis])

        # Outside water, we fix a distance that is bigger than any other
        # distance
        distances[~cut_mask] = max_radius * max_radius + 1

        return np.unravel_index(
            np.argmin(distances, axis=None), distances.shape
        )

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
        return self[jk_m + 1, :, :]

    def bathymetry_in_cells(self) -> np.ndarray:
        """
        Returns a 2d map that for each columns associates the number of water
        cells that are present on that column.

        Returns:
            A 2d numpy array of integers
        """
        return np.count_nonzero(self._mask_array, axis=0)

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
        Best evalutation, since it takes in account e3t.

        Returns:
          a 2d numpy array of floats
        """
        bathymetry = np.sum(self.e3t, axis=0, where=self)
        bathymetry[~self[0, :, :]] = 1e20
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
            file_pointer.variables[mask_var_name], dtype=bool
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


class RegularMask(Mask, RegularMesh):
    def __init__(
        self,
        grid: Grid,
        zlevels: np.ndarray,
        mask_array: np.ndarray,
        e3t: Optional[np.ndarray] = None,
    ):
        if not grid.is_regular():
            raise ValueError("Grid must be regular")

        super().__init__(
            grid=grid, zlevels=zlevels, mask_array=mask_array, e3t=e3t
        )


class Mask(object):
    """
    Defines a mask from a NetCDF file
    """

    def __init__(
        self,
        filename,
        maskvarname="tmask",
        zlevelsvar="nav_lev",
        ylevelsmatvar="nav_lat",
        xlevelsmatvar="nav_lon",
        dzvarname="e3t",
        loadtmask=True,
    ):
        filename = str(filename)

        with netCDF4.Dataset(filename, "r") as dset:
            if (maskvarname in dset.variables) and loadtmask:
                m = dset.variables[maskvarname]
                if len(m.shape) == 4:
                    self._mask = np.array(m[0, :, :, :], dtype=bool)
                elif len(m.shape) == 3:
                    self._mask = np.array(m[:, :, :], dtype=bool)
                else:
                    raise ValueError("Wrong shape: %s" % (m.shape,))
                self._shape = self._mask.shape
            else:
                if loadtmask:
                    raise ValueError(
                        "maskvarname '%s' not found" % (str(maskvarname),)
                    )
                else:
                    dims = dset.dimensions
                    self._shape = (
                        dims["z"].size,
                        dims["y"].size,
                        dims["x"].size,
                    )
            if zlevelsvar in dset.variables:
                if zlevelsvar == "nav_lev":
                    z = dset.variables[zlevelsvar]
                else:
                    z = dset.variables[zlevelsvar][0, :, 0, 0]
                self._zlevels = np.array(z)

                if len(z.shape) != 1:
                    raise ValueError("zlevelsvar must have only one dimension")
                if z.shape[0] not in self._shape:
                    raise ValueError(
                        "cannot match %s lenght with any of %s dimensions"
                        % (zlevelsvar, maskvarname)
                    )
            else:
                raise ValueError(
                    "zlevelsvar '%s' not found" % (str(zlevelsvar),)
                )
            if dzvarname in dset.variables:
                self.e3t = np.array(dset.variables[dzvarname][0, :, :, :])
                # self._dz = np.array(dset.variables[dzvarname][0,:,0,0])
            else:
                if "e3t_0" in dset.variables:
                    self.e3t = np.array(dset.variables["e3t_0"][0, :, :, :])
                else:
                    raise ValueError(
                        "dzvarname '%s' not found" % (str(dzvarname),)
                    )
            if ylevelsmatvar in dset.variables:
                if ylevelsmatvar == "nav_lat":
                    self._ylevels = np.array(dset.variables[ylevelsmatvar])
                if ylevelsmatvar in ["gphit", "gphif", "gphiu", "gphiv"]:
                    self._ylevels = np.array(
                        dset.variables[ylevelsmatvar][0, 0, :, :]
                    )
            else:
                raise ValueError(
                    "ylevelsmatvar '%s' not found" % (str(ylevelsmatvar),)
                )
            if xlevelsmatvar in dset.variables:
                if xlevelsmatvar == "nav_lon":
                    self._xlevels = np.array(dset.variables[xlevelsmatvar])
                if xlevelsmatvar in ["glamt", "glamf", "glamu", "glamv"]:
                    self._xlevels = np.array(
                        dset.variables[xlevelsmatvar][0, 0, :, :]
                    )
            else:
                raise ValueError(
                    "xlevelsmatvar '%s' not found" % (str(xlevelsmatvar),)
                )
            m = dset.variables["e1t"]
            if len(m.shape) == 4:
                self.e1t = np.array(dset.variables["e1t"][0, 0, :, :]).astype(
                    np.float32
                )
                self.e2t = np.array(dset.variables["e2t"][0, 0, :, :]).astype(
                    np.float32
                )
            else:
                self.e1t = np.array(dset.variables["e1t"][0, :, :]).astype(
                    np.float32
                )
                self.e2t = np.array(dset.variables["e2t"][0, :, :]).astype(
                    np.float32
                )

            if loadtmask:
                k = 1
                # Here we look for the first layer from the bottom (i.e. the
                # deepest one) that contains a water cell. In this way, we
                # are sure that in the following part of the code we will read
                # e3t on the deepest column that we have
                while np.all(np.logical_not(self._mask[-k, :])):
                    k += 1
                water_points = np.where(self._mask[-k, :])
                water_point_x, water_point_y = (k[0] for k in water_points)
            else:
                k = 1
                water_point_x, water_point_y = 0, 0

            self._area = self.e1t * self.e2t
            self._dz = self.e3t[:, water_point_x, water_point_y]

            # If we have some layers that contain only land, we use the height
            # of the last land for them
            if k != 1:
                self._dz[-k + 1 :] = self._dz[-k]

        self._regular = None

    def bathymetry(self):
        """
        Calculates the bathymetry used by the model
        Best evalutation, since it takes in account e3t.

        Returns:
        * bathy * a 2d numpy array of floats
        """
        if self.e3t.shape != self.shape:
            print(
                "Warning: e3t is not provided as 3D array in maskfile: Bathymetry will be calculated as function of tmask and zlevels "
            )
            return self.rough_bathymetry()

        cells_bathy = self.bathymetry_in_cells()
        _, jpj, jpi = self.shape
        Bathy = np.ones((jpj, jpi), np.float32) * 1.0e20
        for ji in range(jpi):
            for jj in range(jpj):
                max_lev = cells_bathy[jj, ji]
                if max_lev > 0:
                    Bathy[jj, ji] = self.e3t[:max_lev, jj, ji].sum()
        return Bathy

    def cut_at_level(self, index):
        """
        Arguments:
        * index * integer, depth index

        Returns copy of the mask object, reduced on a single horizontal slice,
        the one of the provided depth index.
        """
        import copy

        New_mask = copy.copy(self)

        _, jpj, jpi = self.shape
        red_mask = np.zeros((1, jpj, jpi), dtype=bool)
        red_mask[0, :, :] = self._mask[index, :, :]
        New_mask._mask = red_mask

        New_mask._shape = red_mask.shape
        New_mask._zlevels = [self._zlevels[index]]
        New_mask._dz = [self._dz[index]]
        return New_mask


class MaskBathymetry(Bathymetry, ABC):
    """
    This class is a bathymetry,  generated starting from a mask, i.e., it
    returns the z-coordinate of the bottom face of the deepest cell of the
    column that contains the point (lon, lat).
    """

    def __init__(self, mask):
        self._mask = mask
        self._bathymetry_data = mask.bathymetry()

        # Fix the bathymetry of the land cells to 0 (to be coherent with the
        # behaviour of the bathymetry classes). Otherwise, if we let the land
        # points to have bathymetry = 1e20, they will be in every depth basin
        self._bathymetry_data[np.logical_not(self._mask.mask[0, :, :])] = 0

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, repr(self._mask))

    def is_inside_domain(self, lon, lat):
        r = 1
        min_lon = self._mask.xlevels.min() - r
        max_lon = self._mask.xlevels.max() + r
        min_lat = self._mask.ylevels.min() - r
        max_lat = self._mask.ylevels.max() + r

        inside_lon = np.logical_and(lon >= min_lon, lon <= max_lon)
        inside_lat = np.logical_and(lat >= min_lat, lat <= max_lat)
        return np.logical_and(inside_lon, inside_lat)

    @staticmethod
    def build_bathymetry(mask):
        if mask.is_regular():
            return RegularMaskBathymetry(mask)
        return NonRegularMaskBathymetry(mask)


class RegularMaskBathymetry(MaskBathymetry):
    def __init__(self, mask):
        if not mask.is_regular():
            raise ValueError(
                "A RegularMaskBathymetry can be generated only from a "
                "regular mask"
            )
        super().__init__(mask)

        assert np.all(
            self._mask.lon[1:] - self._mask.lon[:-1] > 0
        ), "lon array is not ordered"
        assert np.all(
            self._mask.lat[1:] - self._mask.lat[:-1] > 0
        ), "lat array is not ordered"

    def __call__(self, lon, lat):
        lon_index = search_closest_sorted(self._mask.lon, lon)
        lat_index = search_closest_sorted(self._mask.lat, lat)
        return self._bathymetry_data[lat_index, lon_index]


class NonRegularMaskBathymetry(MaskBathymetry):
    def __init__(self, mask):
        super().__init__(mask)
        data = np.stack((self._mask.xlevels, self._mask.ylevels), axis=-1)
        data = data.reshape(-1, 2)

        self._kdtree = spatial.KDTree(data)

    def __call__(self, lon, lat):
        return_float = True
        if hasattr(lon, "__len__") or hasattr(lat, "__len__"):
            return_float = False

        lon, lat = np.broadcast_arrays(lon, lat)

        query_data = np.stack((lon, lat), axis=-1)
        assert query_data.shape[-1] == 2

        _, center_indices = self._kdtree.query(query_data)

        if return_float:
            assert center_indices.shape == (1,) or center_indices.shape == ()
            if center_indices.shape == (1,):
                center_indices = center_indices[0]

        return self._bathymetry_data.flatten()[center_indices]


if __name__ == "__main__":
    Lon = np.arange(10, 50, 0.5)
    Lat = np.arange(30, 46, 0.25)
    L = Grid(Lon, Lat)
    L.convert_lon_lat_to_indices(31.25, 42.02)
    import sys

    sys.exit()
    # Test of convert_lon_lat_wetpoint_indices
    filename = (
        "/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/Somot/meshmask_843_S.nc"
    )
    TheMask = Mask("/g100_work/OGS21_PRACE_P/CLIMA_100/meshmask.nc")
    lat = 33.25925
    lon = 11.18359
    print(TheMask.convert_lon_lat_wetpoint_indices(lon, lat, 2))
    import sys

    sys.exit()
    # TheMask = Mask(filename,zlevelsvar='gdepw', xlevelsmatvar='glamf')
    print(TheMask.is_regular())

    lon = 18.1398
    lat = 37.9585
    i, j = TheMask.convert_lon_lat_wetpoint_indices(lon, lat, 2)
    id, jd = TheMask.convert_lon_lat_wetpoint_indices(lon, lat)
    ip, jp = TheMask.convert_lon_lat_to_indices(lon, lat)
    lon = 9.44
    lat = 40.25
    il, jl = TheMask.convert_lon_lat_wetpoint_indices(lon, lat, 30)
    it, jt = TheMask.convert_lon_lat_wetpoint_indices(lon, lat)
    ipt, jpt = TheMask.convert_lon_lat_to_indices(lon, lat)
    x, y = TheMask.coastline(200, min_line_length=20)
