from commons.helpers import is_number

class Segment(object):
    """Holds a segment.
    """
    def __init__(self, lon_lat_min, lon_lat_max):
        """Segment constructor.

        Args:
            - *lon_lat_min*: tuple or list with Segment's starting point in
              longitude and latitude.
            - *lon_lat_max*: tuple or list with Segment's ending point in
              longitude and latitude.
        """
        if not (isinstance(lon_lat_min, (list, tuple)) and isinstance(lon_lat_max, (list, tuple))):
            raise ValueError("lon_lat_min and lon_lat_max must be tuples or lists")
        if not ((len(lon_lat_min) == 2) and (len(lon_lat_max) == 2)):
            raise ValueError("lon_lat_min and lon_lat_max must hold two values")
        for t in [lon_lat_min, lon_lat_max]:
            for el in t:
                if not is_number(el):
                    raise ValueError("%s is not a number" % (el,))
        self.__lon_min = lon_lat_min[0]
        self.__lat_min = lon_lat_min[1]
        self.__lon_max = lon_lat_max[0]
        self.__lat_max = lon_lat_max[1]

    @property
    def lon_min(self):
        return self.__lon_min

    @property
    def lon_max(self):
        return self.__lon_max

    @property
    def lat_min(self):
        return self.__lat_min

    @property
    def lat_max(self):
        return self.__lat_max

    @property
    def lon_lat_min(self):
        return (self.__lon_min, self.__lat_min)

    @property
    def lon_lat_max(self):
        return (self.__lon_max, self.__lat_max)