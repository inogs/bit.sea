# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
# Utility functions for Regions and Basins

from basin import SimpleBasin
from region import Rectangle

def divide_in_sections(rect, lon_side, lat_side, start_lon, start_lat):
    """
    Cuts a Rectangle Region in Rectangle sections.
    The cutting starts from a point and proceeds to cut along longitude then
    along latitude.

    ++++++++    ++++++++    s->+++++
    ++++++++ -> s->+++++ -> ssssssss
    s->+++++    ssssssss    ssssssss

    Args:
        - *rect*: the Rectangle object to divide.
        - *lon_side*: length in degrees of longitude of a single section.
        - *lat_side*: length in degrees of latitude of a single section.
        - *start_lon*: starting longitude.
        - *start_lat*: starting latitude.

    Returns: a list of SimpleBasin objects.
    """
    # Input validation
    if not isinstance(rect, Rectangle):
        raise ValueError("rect must be an instance of Rectangle")
    if (lon_side > (rect.lonmax - rect.lonmin)) or (lon_side <= 0):
        raise ValueError("Invalid lon_side: %g" % lon_side)
    if (lat_side > (rect.latmax - rect.latmin)) or (lat_side <= 0):
        raise ValueError("Invalid lat_side: %g" % lat_side)
    if not rect.is_inside(start_lon, start_lat):
        raise ValueError("Invalid starting point: %g %g" % (start_lon, start_lat))
    output = list()
    # Bottom Left point
    BL_point = [start_lon, start_lat]
    # Top Right point
    TR_point = [start_lon + lon_side, start_lat + lat_side]
    # Section indices
    lon_in = 0
    lat_in = 0
    # Latitude cycle
    while (TR_point[1] <= rect.latmax):
        # Longitude cycle
        while (TR_point[0] <= rect.lonmax):
            # Create the section
            s = Rectangle(BL_point[0], TR_point[0], BL_point[1], TR_point[1])
            # Create and append the basin to output
            output.append(SimpleBasin("section_%d_%d" % (lat_in, lon_in), rect))
            # Increment longitude index
            lon_in += 1
            # Increment longitude
            BL_point[0] += lon_side
            TR_point[0] += lon_side
        # Reset longitude
        BL_point[0] = start_lon
        TR_point[0] = start_lon + lon_side
        # Reset longitude index
        lon_in = 0
        # Increment latitude index
        lat_in += 1
        # Increment latitude
        BL_point[1] += lat_side
        TR_point[1] += lat_side
    return output
