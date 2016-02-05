# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from commons.utils import is_number

class Segment(object):
    """Holds a segment.
    """
    def __init__(self, lon_lat_min, lon_lat_max, name=None, points=50):
        """Segment constructor.

        Args:
            - *lon_lat_min*: tuple or list with Segment's starting point in
              longitude and latitude.
            - *lon_lat_max*: tuple or list with Segment's ending point in
              longitude and latitude.
            - *name* (optional): a string defining the name of the segment (default: None).
            - *points* (optional): the number of points of the segment
              (default: 50). It is a required argument for diagonals and it is
              ignored for fixed latitude and fixed longitude segments.
        """
        if name is None:
            self.__name = None
        else:
            self.__name = str(name)
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
        self.points = None
        #If the segment is diagonal check that points is a number greater than 0
        if (self.__lon_min != self.__lon_max) and (self.__lat_min != self.__lat_max):
            if not is_number(points) or (points <= 0):
                raise ValueError("points must be defined as a number greater than zero")
            else:
                self.points = points

    def __str__(self):
        if self.__name is None:
            return "Segment from %g,%g to %g,%g" % (self.__lon_min, self.__lat_min, self.__lon_max, self.__lat_max)
        else:
            return "Segment %s from %g,%g to %g,%g" % (self.__name, self.__lon_min, self.__lat_min, self.__lon_max, self.__lat_max)

    def __repr__(self):
        if self.__name is None:
            return "Segment((%f,%f), (%f,%f))" % (self.__lon_min, self.__lat_min, self.__lon_max, self.__lat_max)
        else:
            return "Segment %s ((%f,%f), (%f,%f))" % (self.__name, self.__lon_min, self.__lat_min, self.__lon_max, self.__lat_max)

    @property
    def name(self):
        return self.__name

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
