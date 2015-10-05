import netCDF4

class Mask(object):
    """
    Defines a mask from a NetCDF file
    """
    def __init__(self, filename, maskvarname="tmask", zlevelsvar="nav_lev"):
        filename = str(filename)
        try:
            dset = netCDF4.Dataset(filename)
            if maskvarname in dset.variables:
                #TODO: should it be converted to a NumPy array?
                self.__mask = dset.variables[maskvarname]
            if zlevelsvar in dset.variables:
                #TODO: should it be converted to a NumPy array?
                self.__zlevels = dset.variables[zlevelsvar]
        except:
            raise

    @property
    def mask(self):
        return self.__mask

    @property
    def zlevels(self):
        return self.__zlevels
