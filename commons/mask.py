# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import netCDF4
import pylab as pl

class Mask(object):
    """
    Defines a mask from a NetCDF file
    """
    def __init__(self, filename, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t", loadtmask=True):
        filename = str(filename)
        try:
            dset = netCDF4.Dataset(filename)
            if (maskvarname in dset.variables) and (loadtmask):
                m = dset.variables[maskvarname]
                if len(m.shape) == 4:
                    self._mask = np.array(m[0,:,:,:], dtype=np.bool)
                elif len(m.shape) == 3:
                    self._mask = np.array(m[:,:,:], dtype=np.bool)
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
                z = dset.variables[zlevelsvar]
                if len(z.shape) != 1:
                    raise ValueError("zlevelsvar must have only one dimension")
                if not z.shape[0] in self._shape:
                    raise ValueError("cannot match %s lenght with any of %s dimensions" % (zlevelsvar, maskvarname))
                self._zlevels = np.array(dset.variables[zlevelsvar])
            else:
                raise ValueError("zlevelsvar '%s' not found" % (str(zlevelsvar),))
            if dzvarname in dset.variables:
                self._dz = np.array(dset.variables[dzvarname][0,:,0,0])
            else:
                if 'e3t_0' in dset.variables:
                    self._dz = np.array(dset.variables['e3t_0'][0,:,0,0])
                else:
                    raise ValueError("dzvarname '%s' not found" % (str(dzvarname),))
            if ylevelsmatvar in dset.variables:
                if ylevelsmatvar =='nav_lat':
                    self._ylevels = np.array(dset.variables[ylevelsmatvar])
                if ylevelsmatvar =='gphit':
                    self._ylevels = np.array(dset.variables[ylevelsmatvar][0,0,:,:])
            else:
                raise ValueError("ylevelsmatvar '%s' not found" % (str(ylevelsmatvar),))
            if xlevelsmatvar in dset.variables:
                if xlevelsmatvar=='nav_lon':
                    self._xlevels = np.array(dset.variables[xlevelsmatvar])
                if xlevelsmatvar=='glamt':
                    self._xlevels = np.array(dset.variables[xlevelsmatvar][0,0,:,:])
            else:
                raise ValueError("xlevelsmatvar '%s' not found" % (str(xlevelsmatvar),))
            e1t = np.array(dset.variables['e1t'][0,0,:,:]).astype(np.float32)
            e2t = np.array(dset.variables['e2t'][0,0,:,:]).astype(np.float32)
            self._area = e1t*e2t
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
            raise ValueError("Invalid longitude value: %f (must be between %f and %f)" % (lon, min_lon, max_lon))
        if lat > max_lat or lat < min_lat:
            raise ValueError("Invalid latitude value: %f (must be between %f and %f)" % (lat, min_lat, max_lat))
        #Longitude distances matrix
        d_lon = np.array(self._xlevels - lon)
        d_lon *= d_lon
        #Latitude distances matrix
        d_lat = np.array(self._ylevels - lat)
        d_lat *= d_lat
        #Compute minimum indices
        min_d_lon = d_lon.min()
        min_d_lat = d_lat.min()
        lon_indices = np.where(d_lon == min_d_lon)
        lat_indices = np.where(d_lat == min_d_lat)
        return lon_indices[1][0], lat_indices[0][0]

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
        #Matrixes of indexes of the Mask
        Ilist = np.arange(self.shape[2])
        II = np.tile(Ilist,(self.shape[1],1))
        Jlist = np.arange(self.shape[1]).T
        JJ = np.tile(Jlist,(self.shape[2],1)).T
        IImask = II[self.mask[0,:,:]]
        JJmask = JJ[self.mask[0,:,:]]
        #Find distances from wet points
        distind = (ip-IImask)**2+(jp-JJmask)**2
        #Limit to distance < maxradius
        indd = distind<=maxradius
        ipnarr = IImask[indd]
        jpnarr = JJmask[indd]
        #Assign the first of the nearest wet points 
        if len(ipnarr)>0:
            newip = ipnarr[0]
            newjp = jpnarr[0]
            return newip,newjp
        #If there aren't wet points with distance < maxradius, assign the non-wet point
        else:
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
        '''
        jk_m = self.getDepthIndex(z)
        level_mask = self.mask[jk_m,:,:].copy()
        return level_mask

    def bathymetry_in_cells(self):
        '''
        Returns a 2d array of integers
        '''
        return self._mask.sum(axis=0)

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
        red_mask = np.zeros((1,jpj,jpi),dtype=np.bool)
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

if __name__ == '__main__':
    #Test of convert_lon_lat_wetpoint_indices
    TheMask = Mask('/pico/scratch/userexternal/ateruzzi/DA_COAST_15/wrkdir/MODEL/meshmask.nc')
    lon = 13.4
    lat = 43.6
    i, j = TheMask.convert_lon_lat_wetpoint_indices(lon,lat,2)
    id, jd = TheMask.convert_lon_lat_wetpoint_indices(lon,lat)
    ip, jp = TheMask.convert_lon_lat_to_indices(lon,lat)
    lon = 9.44
    lat = 40.25
    il, jl = TheMask.convert_lon_lat_wetpoint_indices(lon,lat,15)
    it, jt = TheMask.convert_lon_lat_wetpoint_indices(lon,lat)
    ipt, jpt = TheMask.convert_lon_lat_to_indices(lon,lat)
    x,y = TheMask.coastline(200, min_line_length=20)



