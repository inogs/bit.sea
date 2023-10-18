# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from __future__ import print_function

import numpy as np
import netCDF4

from commons.bathymetry import Bathymetry


class OutsideMaskDomain(ValueError):
    pass


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
            self._area = self.e1t*self.e2t
            self._dz = self.e3t[:,0,0]
        except:
            raise

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
            raise OutsideMaskDomain(
                "Invalid longitude value: {} (must be between {} and "
                "{})".format(lon, min_lon, max_lon)
            )
        if lat > max_lat or lat < min_lat:
            raise OutsideMaskDomain(
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

    def mask_at_level(self,z):
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
        if z<self.zlevels[0]: return self.mask[0,:,:].copy()

        jk_m = self.getDepthIndex(z)
        level_mask = self.mask[jk_m+1,:,:].copy()
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
        Calculates a mesh-dependent coastline at a choosen depth.

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
        x1d_0 = self._xlevels[0,:]
        y1d_0 = self._ylevels[:,0]
        X2D,Y2D = np.meshgrid(x1d_0, y1d_0)
        dist = ((X2D - self.xlevels)**2 + (Y2D - self.ylevels)**2).sum()
        regular = dist == 0
        return regular


class MaskBathymetry(Bathymetry):
    """
    This class is a bathymetry,  generated starting from a mask, i.e., it
    returns the z-coordinate of the bottom face of the deepest cell of the
    column that contains the point (lon, lat).
    """
    def __init__(self, mask):
        self._mask = mask
        self._bathymetry_data = mask.bathymetry()

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

    def __call__(self, lon, lat):
        # This is the case when lon and lat are just numbers
        if not hasattr(lon, "__len__") and not hasattr(lat, "__len__"):
            lon_index, lat_index = \
                self._mask.convert_lon_lat_to_indices(lon, lat)
            return self._bathymetry_data[lon_index, lat_index]

        # Now the more complicated case: a bathymetry must be able to handle
        # also numpy arrays. Unfortunately, the method we have used before,
        # convert_lon_lat_to_indices, accepts only two floats.
        # Therefore, we need to perform an iteration. This is not as trivial as
        # it may seem, because we need to take into account the fact that the
        # numpy array may broadcast or the fact that one element may be a numpy
        # array but the other could be a number. First of all, we need to
        # find the right dtype for the output vector.
        data_dtypes = []
        if hasattr(lon, "__len__"):
            data_dtypes.append(np.asarray(lon).dtype)
        else:
            data_dtypes.append(lon)

        if hasattr(lat, "__len__"):
            data_dtypes.append(np.asarray(lat).dtype)
        else:
            data_dtypes.append(lat)
        result_dtype = np.result_type(*data_dtypes)

        # Now we check the right shape
        broadcast_shape = np.broadcast(lon, lat).shape

        lon_broadcast = np.empty(shape=broadcast_shape, dtype=result_dtype)
        lon_broadcast[:] = lon

        lat_broadcast = np.empty(shape=broadcast_shape, dtype=result_dtype)
        lat_broadcast[:] = lat

        # This is the vector that we will return as output
        output = np.empty(
            shape=broadcast_shape,
            dtype=self._bathymetry_data.dtype
        )

        # And this is an iterator over the indices of our arrays
        shape_iterator = np.nditer(output, flags=['multi_index'])
        for _ in shape_iterator:
            current_position = shape_iterator.multi_index

            current_lon = lon_broadcast[current_position]
            current_lat = lat_broadcast[current_position]

            try:
                lon_index, lat_index = self._mask.convert_lon_lat_to_indices(
                    current_lon,
                    current_lat
                )

                output[current_position] = \
                    self._bathymetry_data[lat_index, lon_index]
            except OutsideMaskDomain:
                output[current_position] = 0.

        return output


if __name__ == '__main__':
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



