# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from __future__ import print_function

from abc import ABC
import numpy as np
from pathlib import Path

from scipy import spatial
from sklearn.neighbors import BallTree
import netCDF4
from typing import Optional
from typing import Tuple
from typing import Union

from bitsea.commons.bathymetry import Bathymetry
from bitsea.commons.utils import search_closest_sorted


class OutsideDomain(ValueError):
    pass


class Grid:
    """
    Represents a 2D discretization of a geographical region, storing the
    positions of cell centers in terms of longitude and latitude.

    Args:
        xlevels (2D array): Longitudes of the cell centers.
        ylevels (2D array): Latitudes of the cell centers.
    """
    def __init__(self, *, xlevels: np.ndarray, ylevels: np.ndarray):
        xlevels = np.asarray(xlevels)
        ylevels = np.asarray(ylevels)

        if len(xlevels.shape) != 2:
            raise ValueError(
                f'The lon array is expected to be a 2D array; its current '
                f'shape is {xlevels.shape}'
            )

        if xlevels.shape != ylevels.shape:
            raise ValueError(
                f'The longitude and latitude arrays are expected to have the '
                f'same shape: currently lon = {xlevels.shape} and lat = '
                f'{ylevels.shape}'
            )

        if xlevels.dtype != ylevels.dtype:
            raise ValueError(
                f'The longitude and latitude arrays are expected to have the '
                f'same dtype: currently lon = {xlevels.dtype} and lat = '
                f'{ylevels.dtype}'
            )

        self._xlevels = xlevels
        self._ylevels = ylevels

        self._xlevels.setflags(write=False)
        self._ylevels.setflags(write=False)

        self._lon_min = np.min(self._xlevels)
        self._lat_min = np.min(self._ylevels)

        self._lon_max = np.max(self._xlevels)
        self._lat_max = np.max(self._ylevels)

        self._average_lon_dist = max(
            np.average(np.abs(self._xlevels[:, 1:] - self._xlevels[:, :-1])),
            np.average(np.abs(self._xlevels[1:, :] - self._xlevels[:-1, :]))
        )
        self._average_lat_dist = max(
            np.average(np.abs(self._ylevels[:, 1:] - self._ylevels[:, :-1])),
            np.average(np.abs(self._ylevels[1:, :] - self._ylevels[:-1, :]))
        )

        # BallTree to find the nearest cell center to a given point.
        # Initialized as None and set up only when needed.
        self._balltree: Optional[BallTree] = None

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

    def is_regular(self):
        """
        Returns "True" if the grid is regular, i.e. if all the columns
        of xlevels and all the rows of ylevels are the same
        """
        return False

    def _initialize_balltree(self):
        """
        Initializes the ball tree for this object, enabling efficient
        nearest-neighbor searches.
        """
        # Here we put ylevels and then xlevels because the balltree
        # needs the first coordinate to be the latitude
        data = np.stack(
            (self.ylevels, self.xlevels),
            axis=-1
        )
        data = data.reshape(-1, 2)

        self._balltree = BallTree(data, metric="haversine")

    def _check_if_inside_rectangle(self,
                                   lon: Union[float, np.ndarray],
                                   lat: Union[float, np.ndarray],
                                   tolerance: Optional[float] = None,
                                   raise_if_outside: bool = True)\
            -> Tuple[np.ndarray, np.ndarray]:
        """
        Checks whether one or more points, defined by their
        coordinates, are within a rectangle that bounds all the cell
        centers of this grid.

        Args:
            lon (float or array): Longitude of the point, or an array
              of longitudes if multiple points are being checked.
            lat (float or array): Latitude of the point, or an array
              of latitudes if multiple points are being checked.
            tolerance (float): Expands the boundary of the rectangle by
              this value.
            raise_if_outside (bool): If True, raises an `OutsideDomain`
              exception if any point is outside the rectangle.

        Returns:
            tuple: Two boolean arrays. The first array is `True` for
            points within the longitude range, and the second is `True`
            for points within the latitude range. Performing a logical
            and on these two arrays we get the points that are inside
            the rectangle
        """

        if tolerance is None:
            r = 2 * max(self._average_lat_dist, self._average_lon_dist)
        else:
            r = tolerance
        min_lon = self._lon_min - r
        max_lon = self._lon_max + r
        min_lat = self._lat_min - r
        max_lat = self._lat_max + r

        outside_lon = np.logical_or(lon > max_lon, lon < min_lon)
        outside_lat = np.logical_or(lat > max_lat, lat < min_lat)

        # Here we perform the broadcast just to be sure that we return numpy
        # arrays (even if the input were two floats) and to be sure that we
        # return two arrays of the same dimension
        outside_lon, outside_lat = np.broadcast_arrays(
            outside_lon,
            outside_lat
        )

        if raise_if_outside:
            if np.any(outside_lon):
                if len(outside_lon.shape) > 0:
                    outside_indices = np.nonzero(outside_lon)
                    first_point = tuple(p[0] for p in outside_indices)
                    raise OutsideDomain(
                        f'Point {first_point} with longitude '
                        f'{lon[first_point]} is outside the domain '
                        f'(domain longitude goes from {self._lon_min} to '
                        f'{self._lon_max})'
                    )
                raise OutsideDomain(
                    f'Point with longitude {lon} (and lat = {lat}) is outside '
                    f'the domain (domain latitude goes from {self._lat_min} '
                    f'to {self._lat_max})'
                )

            if np.any(outside_lat):
                if len(outside_lat.shape) > 0:
                    outside_indices = np.nonzero(outside_lat)
                    first_point = tuple(p[0] for p in outside_indices)
                    raise OutsideDomain(
                        f'Point {first_point} with longitude '
                        f'{lon[first_point]} is outside the domain '
                        f'(domain latitude goes from {self._lat_min} to '
                        f'{self._lat_max})'
                    )
                raise OutsideDomain(
                    f'Point with longitude {lon} (and lat = {lat}) is outside '
                    f'the domain (domain latitude goes from {self._lat_min} '
                    f'to {self._lat_max})'
                )

        return outside_lon, outside_lat

    def convert_lon_lat_to_indices(self, *,
                                   lon:Union[float, np.ndarray],
                                   lat:Union[float, np.ndarray]) -> Tuple:
        """Converts longitude and latitude to the nearest indices on
        the mask.

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
        Returns:
            a tuple with two numbers i and j, so that `xlevels[i, j]`
            is the longitude of the closest center point of the grid
            and `ylevels[i, j]` is the latitude of the same point.
        """
        # If lon and lat are numpy array, check if their shape is compatible
        # by simply trying a broadcast
        np.broadcast(lon, lat)

        # Raise an OutsideDomain error if a point is too far from the grid
        self._check_if_inside_rectangle(lon, lat)

        lon, lat = np.broadcast_arrays(lon, lat)

        # Ensure that the balltree is initialized
        if self._balltree is None:
            self._initialize_balltree()

        # Reshape the points so that it becomes a couple of longitudes and
        # latitudes
        query_data = np.stack((lat, lon), axis=-1)
        assert query_data.shape[-1] == 2

        # BallTree wants two-dimensional arrays
        if len(query_data.shape) == 1:
            query_data = query_data[np.newaxis, :]
        if len(query_data.shape) > 2:
            query_data = query_data.reshape((-1, 2))

        _, center_indices = self._balltree.query(query_data)
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


    def convert_i_j_to_lon_lat(self, i: Union[int, np.ndarray],
                               j: Union[int, np.ndarray]) -> Tuple:
        """
        Given two indices, return the corresponding longitude and latitude
        of the point identified by indices i and j.

        Args:
            i (int): The index representing the row.
            j (int): The index representing the column.

        Returns:
            tuple: A tuple containing the longitude and latitude of the point.
        """
        return self._xlevels[i, j], self._ylevels[i, j]

    @staticmethod
    def from_file(file_name: Path, ylevels_var_name: str = "nav_lat",
                  xlevels_var_name: str = "nav_lon"):
        """
        Reads a NetCDF file and returns a `Grid` object.

        Args:
            file_name (Path): Path to the NetCDF file.
            ylevels_var_name (str): Name of the variable in the NetCDF
              file that contains the y-levels (latitude values).
            xlevels_var_name (str): Name of the variable in the NetCDF
              file that contains the x-levels (longitude values).

        Returns:
            Grid: A `Grid` object. If the file contains a regular grid,
            a `RegularGrid` object is returned instead.
        """
        with netCDF4.Dataset(file_name) as f:
            ylevels = np.array(f.variables[ylevels_var_name])
            xlevels = np.array(f.variables[xlevels_var_name])

            if len(xlevels.shape) == 4:
                xlevels = xlevels[0, 0, :, :]
            if len(ylevels.shape) == 4:
                ylevels = ylevels[0, 0, :, :]

        # Now we check if the grid is regular by checking if replicating
        # the first row of the longitudes and the first column of the
        # latitudes we get the same grid that we have read from the file
        x1d = xlevels[0, :]
        y1d = ylevels[:, 0]
        x2d = x1d[np.newaxis, :]
        y2d = y1d[:, np.newaxis]
        dist = np.max((x2d - xlevels) ** 2 + (y2d - ylevels) ** 2)
        if dist < 1e-8:
            return RegularGrid(lon=x1d, lat=y1d)

        return Grid(xlevels=xlevels, ylevels=ylevels)


class RegularGrid(Grid):
    """
    A `RegularGrid` is defined by two vectors: `lon` (longitude) and
    `lat` (latitude).
    It behaves like a standard `Grid`, where the grid points are
    defined as:

        - `xlevels[i, j] = lon[j]` for every `i`
        - `ylevels[i, j] = lat[i]` for every `j`
    """
    def __init__(self, *, lon, lat):
        lon = np.asarray(lon)
        lat = np.asarray(lat)

        if len(lon.shape) != 1:
            raise ValueError(
                f"lon should be a 1D array: its current shape is {lon.shape}"
            )
        if len(lat.shape) != 1:
            raise ValueError(
                f"lat should be a 1D array: its current shape is {lat.shape}"
            )

        grid_shape = (lat.shape[0], lon.shape[0])
        xlevels = np.broadcast_to(lon, grid_shape)
        ylevels = np.broadcast_to(lat[:, np.newaxis], grid_shape)

        super().__init__(xlevels=xlevels, ylevels=ylevels)

    @property
    def lon(self):
        return self._xlevels[0, :]

    @property
    def lat(self):
        return self._ylevels[:, 0]

    def is_regular(self):
        return True

    def convert_i_j_to_lon_lat(self, i, j):
        return self.lon[j], self.lat[i]

    def convert_lon_lat_to_indices(self, *, lon, lat):
        return_array = False
        if isinstance(lon, np.ndarray) or isinstance(lat, np.ndarray):
            return_array = True

        lon_distances = np.abs(self.lon[:, np.newaxis] - lon)
        lat_distances = np.abs(self.lat[:, np.newaxis] - lat)

        i = np.argmin(lon_distances, axis=0)
        j = np.argmin(lat_distances, axis=0)

        if not return_array:
            i = int(i)
            j = int(j)

        return i, j


class Mask(object):
    """
    Defines a mask from a NetCDF file
    """
    def __init__(self, filename, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t", loadtmask=True):
        filename = str(filename)

        with netCDF4.Dataset(filename, 'r') as dset:
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
                    raise ValueError("maskvarname '%s' not found" % (str(maskvarname),))
                else:
                    dims = dset.dimensions
                    self._shape = (dims['z'].size, dims['y'].size, dims['x'].size)
            if zlevelsvar in dset.variables:
                if zlevelsvar == "nav_lev":
                    z = dset.variables[zlevelsvar]
                else:
                    z = dset.variables[zlevelsvar][0,:,0,0]
                self._zlevels = np.array(z)

                if len(z.shape) != 1:
                    raise ValueError("zlevelsvar must have only one dimension")
                if not z.shape[0] in self._shape:
                    raise ValueError("cannot match %s lenght with any of %s dimensions" % (zlevelsvar, maskvarname))
            else:
                raise ValueError("zlevelsvar '%s' not found" % (str(zlevelsvar),))
            if dzvarname in dset.variables:
                self.e3t =  np.array(dset.variables[dzvarname][0,:,:,:])
                #self._dz = np.array(dset.variables[dzvarname][0,:,0,0])
            else:
                if 'e3t_0' in dset.variables:
                    self.e3t = np.array(dset.variables['e3t_0'][0,:,:,:])
                else:
                    raise ValueError("dzvarname '%s' not found" % (str(dzvarname),))
            if ylevelsmatvar in dset.variables:
                if ylevelsmatvar =='nav_lat': self._ylevels = np.array(dset.variables[ylevelsmatvar])
                if ylevelsmatvar in ['gphit','gphif','gphiu','gphiv']: self._ylevels = np.array(dset.variables[ylevelsmatvar][0,0,:,:])
            else:
                raise ValueError("ylevelsmatvar '%s' not found" % (str(ylevelsmatvar),))
            if xlevelsmatvar in dset.variables:
                if xlevelsmatvar=='nav_lon': self._xlevels = np.array(dset.variables[xlevelsmatvar])
                if xlevelsmatvar in ['glamt','glamf','glamu', 'glamv']: self._xlevels = np.array(dset.variables[xlevelsmatvar][0,0,:,:])
            else:
                raise ValueError("xlevelsmatvar '%s' not found" % (str(xlevelsmatvar),))
            m = dset.variables['e1t']
            if len(m.shape) == 4:
                self.e1t = np.array(dset.variables['e1t'][0,0,:,:]).astype(np.float32)
                self.e2t = np.array(dset.variables['e2t'][0,0,:,:]).astype(np.float32)
            else:
                self.e1t = np.array(dset.variables['e1t'][0,:,:]).astype(np.float32)
                self.e2t = np.array(dset.variables['e2t'][0,:,:]).astype(np.float32)

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
                self._dz[-k + 1:] = self._dz[-k]

        self._regular = None

    @property
    def mask(self):
        return self._mask

    @property
    def xlevels(self):
        return self._xlevels

    @property
    def ylevels(self):
        return self._ylevels

    @property
    def zlevels(self):
        return self._zlevels
    @property
    def lon(self):
        return self._xlevels[0,:]
    @property
    def lat(self):
        return self._ylevels[:,0]
    @property
    def jpi(self):
        return self._shape[2]
    @property
    def jpj(self):
       return self._shape[1]
    @property
    def jpk(self):
        return self._shape[0]
    @property
    def dz(self):
        return self._dz

    @property
    def shape(self):
        return self._shape
    @property
    def area(self):
        return self._area

    def convert_lon_lat_to_indices(self, lon, lat):
        """Converts longitude and latitude to the nearest indices on the mask.

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
        Returns: a tuple of numbers, the first one is the longitude index and
        the other one is the latitude index.
        """
        #Input validation
        lon = float(lon)
        lat = float(lat)
        r=1.0
        min_lon = self._xlevels.min()-r
        max_lon = self._xlevels.max()+r
        min_lat = self._ylevels.min()-r
        max_lat = self._ylevels.max()+r
        if lon > max_lon or lon < min_lon:
            raise OutsideDomain(
                "Invalid longitude value: {} (must be between {} and "
                "{})".format(lon, min_lon, max_lon)
            )
        if lat > max_lat or lat < min_lat:
            raise OutsideDomain(
                "Invalid latitude value: {} (must be between {} and "
                "{})".format(lat, min_lat, max_lat)
            )

        #Longitude distances matrix
        d_lon = np.array(self._xlevels - lon)
        d_lon *= d_lon
        #Latitude distances matrix
        d_lat = np.array(self._ylevels - lat)
        d_lat *= d_lat

        if self.is_regular():
            #Compute minimum indices
            min_d_lon = d_lon.min()
            min_d_lat = d_lat.min()
            lon_index = np.where(d_lon == min_d_lon)[1][0]
            lat_index = np.where(d_lat == min_d_lat)[0][0]
        else:
            dist= d_lon + d_lat
            indlat, indlon =np.where(dist==dist.min())
            lon_index=indlon[0]
            lat_index=indlat[0]

        return lon_index, lat_index



    def convert_lon_lat_wetpoint_indices(self, lon, lat, maxradius=2):
        """Converts longitude and latitude to the nearest water point indices on the mask with maximum distance limit

        Args:
            - *lon*: Longitude in degrees.
            - *lat*: Latitude in degrees.
            - *maxradius* : Maximum distance where the water point is searched (in grid units, integer, default: 2)
        Returns: a tuple of numbers, the first one is the longitude index and
        the other one is the latitude index.
        """
        #Indexes of the input lon, lat
        lon = float(lon)
        lat = float(lat)
        ip,jp = self.convert_lon_lat_to_indices(lon,lat)
        if self.mask[0,jp,ip] : return ip, jp
        #Matrixes of indexes of the Mask
        Ilist = np.arange(self.shape[2])
        II = np.tile(Ilist,(self.shape[1],1))
        Jlist = np.arange(self.shape[1]).T
        JJ = np.tile(Jlist,(self.shape[2],1)).T
        IImask = II[self.mask[0,:,:]]
        JJmask = JJ[self.mask[0,:,:]]
        #Find distances from wet points
        distind = (ip-IImask)**2+(jp-JJmask)**2
        if maxradius is None :
            maxradius = distind.min()
        indd = distind<=maxradius
        #Limit to distance < maxradius)
        ipnarr = IImask[indd]
        jpnarr = JJmask[indd]
        #Assign the nearest wet points 
        if len(ipnarr)>0:
            distmask = distind[indd]
            indmin = np.argmin(distmask)
            newip = ipnarr[indmin]
            newjp = jpnarr[indmin]
            return newip,newjp
        #If there aren't wet points with distance < maxradius, assign the non-wet point
        else:
            print('WARNING: Using terrain point indexes, put maxradius=',distind.min(), ' or maxradius=None')
            return ip,jp


    def convert_i_j_to_lon_lat(self, i, j):
        """Converts i and j indexes to longitude and latitude of center of cells
        i is indeded as longitudinal index, as well as j is latitudinal index
        """
        return (self._xlevels[j,i], self._ylevels[j,i])

    def getDepthIndex(self, z):
        '''Converts a depth expressed in meters in the corresponding index level
        The returned value is an integer indicating the previous (not nearest) depth in the z-levels.

        Example:
        M = Mask(filename)
        k = M.getDepthIndex(200.)
        M.zlevels[k]
          returns 192.60
        '''
        jk_m = 0
        for jk,depth in enumerate(self.zlevels):
            if depth < z:
                jk_m=jk
        return jk_m

    def mask_at_level(self, z):
        '''
        Returns a 2d map of logicals, for the depth (m) provided as argument.
        as a slice of the tmask.

        When z is not a point in the discretization, jk_m is selected as the
        point immediately before. This depth level should not be included in
        the mask:

        (jk_m-1)   FFFFFFFFFFFFFFFF
        (jk_m  )   FFFFFFFFFFFFFFFF
        z----------------------
        (jk_m +1)  TTTTTTTTTTTTTTTT
        (jk_m +2)  TTTTTTTTTTTTTTTT
        '''
        if z < self.zlevels[0]:
            return self.mask[0, :, :].copy()

        jk_m = self.getDepthIndex(z)
        level_mask = self.mask[jk_m+1, :, :].copy()
        return level_mask

    def bathymetry_in_cells(self):
        '''
        Returns a 2d array of integers
        '''
        return self._mask.sum(axis=0)
    def rough_bathymetry(self):
        '''
        Calculates the bathymetry used by the model 
        It does not not takes in account e3t

        Returns:
        * bathy * a 2d numpy array of floats
        '''
        Cells = self.bathymetry_in_cells()
        zlevels =np.concatenate((np.array([0]) , self.zlevels))
        bathy = zlevels[Cells]
        return bathy

    def bathymetry(self):
        '''
        Calculates the bathymetry used by the model
        Best evalutation, since it takes in account e3t.

        Returns:
        * bathy * a 2d numpy array of floats
        '''
        if (self.e3t.shape !=self.shape ) :
            print("Warning: e3t is not provided as 3D array in maskfile: Bathymetry will be calculated as function of tmask and zlevels ")
            return self.rough_bathymetry()

        cells_bathy = self.bathymetry_in_cells()
        _,jpj,jpi = self.shape
        Bathy = np.ones((jpj,jpi),np.float32)*1.e+20
        for ji in range(jpi):
            for jj in range(jpj):
                max_lev=cells_bathy[jj,ji]
                if max_lev > 0:
                    Bathy[jj,ji] = self.e3t[:max_lev,jj,ji].sum()
        return Bathy

    def cut_at_level(self,index):
        '''
        Arguments:
        * index * integer, depth index

        Returns copy of the mask object, reduced on a single horizontal slice,
        the one of the provided depth index.
        '''
        import copy
        New_mask = copy.copy(self)

        _,jpj,jpi = self.shape
        red_mask = np.zeros((1,jpj,jpi),dtype=bool)
        red_mask[0,:,:] = self._mask[index,:,:]
        New_mask._mask = red_mask

        New_mask._shape = red_mask.shape
        New_mask._zlevels = [self._zlevels[index]]
        New_mask._dz      = [self._dz[index]]
        return New_mask
    def coastline(self,depth, min_line_length=30):
        '''
        Calculates a mesh-dependent coastline at a chosen depth.

        Arguments :
        * level           * depths expressed in meters
        * min_line_length * integer indicating the minimum number of points of each coastline,
                            it is a threshold to avoid a lot of small islands.
        Returns:
        * x,y *  numpy 1d arrays, containing nans to separate the lines, in order to be easily plotted.
        '''
        import matplotlib.pyplot as pl 
        tmask= self.mask_at_level(depth).astype(np.float64)

        H = pl.contour(self.xlevels, self.ylevels, tmask, levels=[float(0.5)])

        PATH_LIST = H.collections[0].get_paths()

        X = np.zeros((0,2))
        nan = np.zeros((1,2))*np.nan
        for p in PATH_LIST:
            v = p.vertices
            nPoints, _ = v.shape
            if nPoints > min_line_length:
                X = np.concatenate((X, v, nan),axis=0)
        x = X[:-1,0]
        y = X[:-1,1]
        return x,y
    def is_regular(self):
        '''
        Returns True if a mesh is regular, False if is not.
        Regular means that (xlevels, ylevels) can be obtained by np.meshgrid(xlevels[k,:], ylevels[:,k])
        '''
        if self._regular is None:
            x1d_0 = self._xlevels[0,:]
            y1d_0 = self._ylevels[:,0]
            X2D, Y2D = np.meshgrid(x1d_0, y1d_0)
            dist = ((X2D - self.xlevels)**2 + (Y2D - self.ylevels)**2).sum()
            self._regular = dist == 0
        return self._regular


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
                'A RegularMaskBathymetry can be generated only from a '
                'regular mask'
            )
        super().__init__(mask)

        assert np.all(self._mask.lon[1:] - self._mask.lon[:-1] > 0), \
            "lon array is not ordered"
        assert np.all(self._mask.lat[1:] - self._mask.lat[:-1] > 0), \
            "lat array is not ordered"

    def __call__(self, lon, lat):
        lon_index = search_closest_sorted(self._mask.lon, lon)
        lat_index = search_closest_sorted(self._mask.lat, lat)
        return self._bathymetry_data[lat_index, lon_index]


class NonRegularMaskBathymetry(MaskBathymetry):
    def __init__(self, mask):
        super().__init__(mask)
        data = np.stack(
            (self._mask.xlevels, self._mask.ylevels),
            axis=-1
        )
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


if __name__ == '__main__':
    Lon= np.arange(10,50,0.5)
    Lat = np.arange(30,46,0.25)
    L = Grid(Lon,Lat)
    L.convert_lon_lat_to_indices(31.25, 42.02)
    import sys
    sys.exit()
    #Test of convert_lon_lat_wetpoint_indices
    filename="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/Somot/meshmask_843_S.nc"
    TheMask = Mask('/g100_work/OGS21_PRACE_P/CLIMA_100/meshmask.nc')
    lat=33.25925
    lon=11.18359
    print(TheMask.convert_lon_lat_wetpoint_indices(lon,lat,2))
    import sys
    sys.exit()
    #TheMask = Mask(filename,zlevelsvar='gdepw', xlevelsmatvar='glamf')
    print(TheMask.is_regular())

    lon = 18.1398
    lat = 37.9585
    i, j = TheMask.convert_lon_lat_wetpoint_indices(lon,lat,2)
    id, jd = TheMask.convert_lon_lat_wetpoint_indices(lon,lat)
    ip, jp = TheMask.convert_lon_lat_to_indices(lon,lat)
    lon = 9.44
    lat = 40.25
    il, jl = TheMask.convert_lon_lat_wetpoint_indices(lon,lat,30)
    it, jt = TheMask.convert_lon_lat_wetpoint_indices(lon,lat)
    ipt, jpt = TheMask.convert_lon_lat_to_indices(lon,lat)
    x,y = TheMask.coastline(200, min_line_length=20)



