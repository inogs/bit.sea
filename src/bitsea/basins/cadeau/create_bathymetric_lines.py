from os import path
from pathlib import Path

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cv2 as cv
import matplotlib.pyplot as plt
import numpy as np
import shapely

from bitsea.commons.mask import Mask, MaskBathymetry
from bitsea.commons.bathymetry import GEBCOBathymetry, SmoothBathymetry
from bitsea.basins.region import BathymetricPolygon, Polygon
from bitsea.basins.basin import SimpleBathymetricBasin

from bitsea.basins.cadeau.nad_V0 import \
    BATHYMETRIC_ISOLINES_DIR, OUTSIDE_NORTH, DOMAIN_LIMIT_SOUTH, \
    OUTSIDE_WEST, DEPTHS


CURRENT_FILE_DIR = Path(path.dirname(path.realpath(__file__)))

BATHYMETRIC_FILE_NAME = 'gebco_2023_n46.8201_s41.0522_w11.0925_e20.0354.nc'
BATHYMETRIC_DATA_FILE = CURRENT_FILE_DIR / 'data' / BATHYMETRIC_FILE_NAME

GEBCO_BATHYMETRY = GEBCOBathymetry(BATHYMETRIC_DATA_FILE)
BATHYMETRY = SmoothBathymetry(GEBCO_BATHYMETRY, 0.03, 10)


DEPTHS = (20, 30, 35)

OUTSIDE_NORTH = 46.
OUTSIDE_WEST = 12.
DOMAIN_LIMIT_SOUTH = 43.47

SIMPLIFY = True
TOLERANCE = 0.001

POINTS = 2_000

MAIN_REGION_LAT = [
    DOMAIN_LIMIT_SOUTH, 44.70, OUTSIDE_NORTH, OUTSIDE_NORTH, DOMAIN_LIMIT_SOUTH
]
MAIN_REGION_LON = [
    15.418, 13.91, 13.91, OUTSIDE_WEST, OUTSIDE_WEST
]


CURRENT_FILE_DIR = Path(path.dirname(path.realpath(__file__)))
BATHYMETRIC_ISOLINES_DIR = CURRENT_FILE_DIR / 'bathymetric_isolines'


def build_main_region_shallower_than(depth):
    bathymetric_polygon = BathymetricPolygon(
        MAIN_REGION_LON,
        MAIN_REGION_LAT,
        BATHYMETRY,
        shallower_than=depth
    )
    return SimpleBathymetricBasin(
        'main_zone_{}'.format(depth),
        bathymetric_polygon
    )


def get_bathymetric_isoline(coords, tolerance=0.0001):
    lon = coords[:, 0]
    lat = coords[:, 1]

    # remove all the coordinates that are outside the domain on the west
    outside = lon < OUTSIDE_WEST + tolerance

    # The same, but on the north
    outside = np.logical_or(
        outside,
        lat > OUTSIDE_NORTH - tolerance
    )

    # Now we remove also the points that are below the south part of the domain
    outside = np.logical_or(
        outside,
        lat < DOMAIN_LIMIT_SOUTH + tolerance
    )

    lat = lat[~outside]
    lon = lon[~outside]

    # Ensure that the first part of the line is near Trieste and not near Ancona
    if lat[0] < lat[-1]:
        lat = lat[::-1]
        lon = lon[::-1]

    # We want the last point to be on the south boundary of the domain
    last_point_lat = DOMAIN_LIMIT_SOUTH
    last_point_lon = lon[-1]

    lat = np.append(lat, last_point_lat)
    lon = np.append(lon, last_point_lon)

    return np.stack([lon, lat]).T


def generate_bathymetric_isolines():
    basins = tuple(
        build_main_region_shallower_than(i) for i in DEPTHS
    )

    for (i, basin), depth in zip(enumerate(basins), DEPTHS):
        lon, lat, data = basin.grid(
            lon_window=(12, 16),
            lat_window=(43, 46),
            lon_points=POINTS,
            lat_points=2 * POINTS,
        )

        tableau_poche = data.astype(np.uint8)
        contours, hierarchy = cv.findContours(tableau_poche, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)

        npoints=0
        for ncon in range(len(contours)):
            npoints_test=contours[ncon][:,0,0].shape[0]
            if npoints_test > npoints:
                x=contours[ncon][:,0,0]
                y=contours[ncon][:,0,1]
                npoints = npoints_test

        if SIMPLIFY:
            poly = shapely.geometry.Polygon(tuple((x, y) for x, y in zip(lon[x], lat[y])))
            poly_s = poly.simplify(tolerance=TOLERANCE)
            coords = np.asarray(poly_s.boundary.coords[:])
        else:
            coords = np.stack((lon[x], lat[y]), axis=-1)

        coords = get_bathymetric_isoline(coords)

        np.savetxt(
            BATHYMETRIC_ISOLINES_DIR / 'depth_{}.txt'.format(depth),
            coords
        )



if __name__ == '__main__':
    generate_bathymetric_isolines()

