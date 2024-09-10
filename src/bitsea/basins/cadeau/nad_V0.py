from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from numbers import Number
from os import path

import numpy as np

from bitsea.basins.region import Polygon, Rectangle
from bitsea.basins.basin import SimpleBasin, SimplePolygonalBasin, ComposedBasin


# If you change these parameters, please re-generate the bathymetric
# lines using the script create_bathymetric_lines.py
DEPTHS = (20, 30, 35)
DEPTH_NAMES = ('costal', 'shallow', 'mid', 'deep')

OUTSIDE_NORTH = 46.
OUTSIDE_WEST = 12.
OUTSIDE_EAST = 16.5
DOMAIN_LIMIT_SOUTH = 43.47


@dataclass
class Point:
    lon: Number
    lat: Number


# Zone parameters
TRIESTE_GULF_SOUTH_LIMIT = 45.496

ISTRIAN_MIDDLE_LINE = 45.2

CENTER_POINT = Point(13.285, 44.61)

PREMANTURA_LONGITUDE = 13.91

SOUTH_MIDDLE_LINE = Point(14.6, 43.47)
SOUTH_EASTERN_LINE = Point(15.418, 43.47)

SOUTH_LAT_CUT_1 = 43.965
SOUTH_LAT_CUT_2 = 43.605

ZONE_B_HORIZONTAL_CUT = 45.11

ZONE_B_R_1_START = Point(12.223, 45.500)
ZONE_B_R_1_END =  Point(12.801, 45.301)

ZONE_B_R_2_START = Point(12.223, 45.223)
ZONE_B_R_2_END = Point(13.074, 45.184)


CURRENT_FILE_DIR = Path(path.dirname(path.realpath(__file__)))
BATHYMETRIC_ISOLINES_DIR = CURRENT_FILE_DIR / 'bathymetric_isolines'


class Line:
    def __init__(self, p1: Point, p2: Point):
        self.p1 = p1
        self.p2 = p2

    def get_point_on_latitude(self, lat):
        if np.abs(self.p1.lat - self.p2.lat) < 1e-7:
            raise ValueError(
                'The line is a parallel; it is not possible to find a point '
                'with a specific latitude'
            )

        segment_position = (lat - self.p1.lat) / (self.p2.lat - self.p1.lat)
        lon = self.p1.lon + (self.p2.lon - self.p1.lon) * segment_position
        return Point(lon, lat)

    def get_point_on_longitude(self, lon):
        if np.abs(self.p1.lon - self.p2.lon) < 1e-7:
            raise ValueError(
                'The line is a meridian; it is not possible to find a point '
                'with a specific longitude'
            )

        segment_position = (lon - self.p1.lon) / (self.p2.lon - self.p1.lon)
        lat = self.p1.lat + (self.p2.lat - self.p1.lat) * segment_position
        return Point(lon, lat)


@lru_cache(1)
def generate_bathymetric_regions():
    lines = {}
    for depth in DEPTHS:
        file_path = BATHYMETRIC_ISOLINES_DIR / 'depth_{}.txt'.format(depth)
        current_lines = np.loadtxt(file_path)
        lines[depth] = current_lines

    regions = {}
    if len(lines) > 0:
        line0 = lines[DEPTHS[0]]
        lons = line0[:, 0]
        lats = line0[:, 1]

        start_polygon_lons = np.array([
            OUTSIDE_WEST, lons[0]
        ])

        start_polygon_lats = np.array([
            OUTSIDE_NORTH, OUTSIDE_NORTH
        ])

        end_polygon_lons = np.array([
            OUTSIDE_WEST
        ])
        end_polygon_lats = np.array([
            DOMAIN_LIMIT_SOUTH
        ])

        lons = np.concatenate((
            start_polygon_lons, lons, end_polygon_lons
        ))
        lats = np.concatenate((
            start_polygon_lats, lats, end_polygon_lats
        ))

        current_region = Polygon(lons, lats)
        regions[DEPTH_NAMES[0]] = current_region

    for i, (d1, d2) in enumerate(zip(DEPTHS[:-1], DEPTHS[1:]), 1):
        line1 = lines[d1]
        line2 = lines[d2]

        lons = np.concatenate((
            line2[:, 0],
            line1[::-1, 0],
        ))

        lats = np.concatenate((
            line2[:, 1],
            line1[::-1, 1],
        ))

        current_region = Polygon(lons, lats)
        regions[DEPTH_NAMES[i]] = current_region

    if len(lines) > 0:
        line_last = lines[DEPTHS[-1]]

        lons = line_last[:, 0]
        lats = line_last[:, 1]

        eastern_line = Line(
            Point(PREMANTURA_LONGITUDE, CENTER_POINT.lat),
            SOUTH_EASTERN_LINE
        )

        end_polygon_lons = np.array([
            eastern_line.get_point_on_latitude(lats[-1]).lon, lons[0]
        ])
        end_polygon_lats = np.array([
            lats[-1], lats[0]
        ])
        lons = np.concatenate((lons, end_polygon_lons))
        lats = np.concatenate((lats, end_polygon_lats))

        current_region = Polygon(lons, lats)
        regions[DEPTH_NAMES[-1]] = current_region

    return regions


def plot_bathymetric_regions():
    import matplotlib.pyplot as plt
    import cartopy
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    n_points = 2_000

    plt.figure(dpi=600, figsize=(16, 9))
    axes = plt.axes(projection=ccrs.PlateCarree())
    regions = generate_bathymetric_regions()
    for i, (title, region) in enumerate(regions.items()):
        basin = SimpleBasin(title, region)
        basin.plot(
            lon_window=(OUTSIDE_EAST, OUTSIDE_WEST),
            lat_window=(DOMAIN_LIMIT_SOUTH - 0.05, OUTSIDE_NORTH),
            lon_points=n_points,
            lat_points=2 * n_points,
            color='C{}'.format(i),
            axes=axes
        )

    land = cfeature.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='10m',
        facecolor=(0.5, 0.5, 0.5, 0.8)
    )

    axes.add_feature(land, zorder=2)
    plt.savefig('bathymetric_isolines.pdf')
    plt.close()


@lru_cache(1)
def generate_zones():
    zone_A_region = Rectangle(
        CENTER_POINT.lon, PREMANTURA_LONGITUDE,
        TRIESTE_GULF_SOUTH_LIMIT, OUTSIDE_NORTH
    )

    zone_B_region = Rectangle(
        OUTSIDE_WEST, CENTER_POINT.lon, CENTER_POINT.lat, OUTSIDE_NORTH
    )

    r1 = Line(ZONE_B_R_1_START, ZONE_B_R_1_END)
    r2 = Line(ZONE_B_R_2_START, ZONE_B_R_2_END)

    zone_B1_region = Polygon(
        (
            CENTER_POINT.lon, CENTER_POINT.lon,
            OUTSIDE_WEST, OUTSIDE_WEST, ZONE_B_R_1_END.lon,
            CENTER_POINT.lon
        ),(
            ZONE_B_R_1_END.lat, OUTSIDE_NORTH,
            OUTSIDE_NORTH, r1.get_point_on_longitude(OUTSIDE_WEST).lat,
            ZONE_B_R_1_END.lat, ZONE_B_R_1_END.lat
        )
    )

    zone_B2_region = Polygon(
        (
            CENTER_POINT.lon, CENTER_POINT.lon, OUTSIDE_WEST,
            OUTSIDE_WEST, CENTER_POINT.lon
        ),(
            ZONE_B_HORIZONTAL_CUT, r2.get_point_on_longitude(CENTER_POINT.lon).lat,
            r2.get_point_on_longitude(OUTSIDE_WEST).lat, ZONE_B_HORIZONTAL_CUT,
            ZONE_B_HORIZONTAL_CUT
        )
    )

    zone_B3_region = Polygon(
        (
            CENTER_POINT.lon, CENTER_POINT.lon,
            ZONE_B_R_1_END.lon, OUTSIDE_WEST,
            OUTSIDE_WEST, CENTER_POINT.lon
        ),(
            r2.get_point_on_longitude(CENTER_POINT.lon).lat,
            ZONE_B_R_1_END.lat, ZONE_B_R_1_END.lat,
            r1.get_point_on_longitude(OUTSIDE_WEST).lat,
            r2.get_point_on_longitude(OUTSIDE_WEST).lat,
            r2.get_point_on_longitude(CENTER_POINT.lon).lat,
        )
    )

    zone_B4_region = Rectangle(
        OUTSIDE_WEST, CENTER_POINT.lon, CENTER_POINT.lat, ZONE_B_HORIZONTAL_CUT
    )

    zone_C_region = Rectangle(
        CENTER_POINT.lon, PREMANTURA_LONGITUDE,
        ISTRIAN_MIDDLE_LINE, TRIESTE_GULF_SOUTH_LIMIT,
    )

    zone_D_region = Rectangle(
        CENTER_POINT.lon, PREMANTURA_LONGITUDE,
        CENTER_POINT.lat, ISTRIAN_MIDDLE_LINE,
    )

    middle_line = Line(CENTER_POINT, SOUTH_MIDDLE_LINE)
    eastern_line = Line(
        Point(PREMANTURA_LONGITUDE, CENTER_POINT.lat),
        SOUTH_EASTERN_LINE
    )

    middle_p_cut_1 = middle_line.get_point_on_latitude(SOUTH_LAT_CUT_1)
    east_p_cut_1 = eastern_line.get_point_on_latitude(SOUTH_LAT_CUT_1)

    zone_E_region = Polygon(
        (
            OUTSIDE_EAST, OUTSIDE_EAST, PREMANTURA_LONGITUDE,
            PREMANTURA_LONGITUDE, east_p_cut_1.lon, OUTSIDE_EAST
        ), (
            east_p_cut_1.lat, OUTSIDE_NORTH, OUTSIDE_NORTH,
            CENTER_POINT.lat, east_p_cut_1.lat, east_p_cut_1.lat
        )
    )
    zone_F_region = Polygon(
        (middle_p_cut_1.lon, CENTER_POINT.lon, OUTSIDE_WEST, OUTSIDE_WEST, middle_p_cut_1.lon),
        (middle_p_cut_1.lat, CENTER_POINT.lat, CENTER_POINT.lat, middle_p_cut_1.lat, middle_p_cut_1.lat)
    )

    zone_G_region = Polygon(
        (east_p_cut_1.lon, PREMANTURA_LONGITUDE, CENTER_POINT.lon, middle_p_cut_1.lon, east_p_cut_1.lon),
        (east_p_cut_1.lat, CENTER_POINT.lat, CENTER_POINT.lat, middle_p_cut_1.lat, east_p_cut_1.lat)
    )

    middle_p_cut_2 = middle_line.get_point_on_latitude(SOUTH_LAT_CUT_2)
    east_p_cut_2 = eastern_line.get_point_on_latitude(SOUTH_LAT_CUT_2)

    zone_H_region = Polygon(
        (middle_p_cut_2.lon, middle_p_cut_1.lon, OUTSIDE_WEST, OUTSIDE_WEST, middle_p_cut_2.lon),
        (middle_p_cut_2.lat, middle_p_cut_1.lat, middle_p_cut_1.lat, middle_p_cut_2.lat, middle_p_cut_2.lat)
    )

    zone_I_region = Polygon(
        (east_p_cut_2.lon, east_p_cut_1.lon, middle_p_cut_1.lon, middle_p_cut_2.lon, east_p_cut_2.lon),
        (east_p_cut_2.lat, east_p_cut_1.lat, middle_p_cut_1.lat, middle_p_cut_2.lat, east_p_cut_2.lat)
    )

    zone_J_region = Polygon(
        (OUTSIDE_EAST, OUTSIDE_EAST, east_p_cut_1.lon, east_p_cut_2.lon, OUTSIDE_EAST),
        (east_p_cut_2.lat, east_p_cut_1.lat, east_p_cut_1.lat, east_p_cut_2.lat, east_p_cut_2.lat)
    )

    middle_p_cut_3 = middle_line.get_point_on_latitude(DOMAIN_LIMIT_SOUTH)
    east_p_cut_3 = eastern_line.get_point_on_latitude(DOMAIN_LIMIT_SOUTH)

    zone_K_region = Polygon(
        (middle_p_cut_3.lon, middle_p_cut_2.lon, OUTSIDE_WEST, OUTSIDE_WEST, middle_p_cut_3.lon),
        (DOMAIN_LIMIT_SOUTH, middle_p_cut_2.lat, middle_p_cut_2.lat, DOMAIN_LIMIT_SOUTH, DOMAIN_LIMIT_SOUTH)
    )

    zone_L_region = Polygon(
        (east_p_cut_3.lon, east_p_cut_2.lon, middle_p_cut_2.lon, middle_p_cut_3.lon, east_p_cut_3.lon),
        (DOMAIN_LIMIT_SOUTH, middle_p_cut_2.lat, middle_p_cut_2.lat, DOMAIN_LIMIT_SOUTH, DOMAIN_LIMIT_SOUTH)
    )

    zone_M_region = Polygon(
        (OUTSIDE_EAST, OUTSIDE_EAST, east_p_cut_2.lon, east_p_cut_3.lon, OUTSIDE_EAST),
        (DOMAIN_LIMIT_SOUTH, middle_p_cut_2.lat, middle_p_cut_2.lat, DOMAIN_LIMIT_SOUTH, DOMAIN_LIMIT_SOUTH)
    )

    zones = (
        SimplePolygonalBasin('zoneA', zone_A_region, 'Gulf of Trieste'),
        SimplePolygonalBasin('zoneB', zone_B_region, 'West North Adri'),
        SimplePolygonalBasin('zoneB1', zone_B1_region, 'Zone B1'),
        SimplePolygonalBasin('zoneB2', zone_B2_region, 'Venice'),
        SimplePolygonalBasin('zoneB3', zone_B3_region, 'Between Venice and Po Delta'),
        SimplePolygonalBasin('zoneB4', zone_B4_region, 'Po ROFI'),
        SimplePolygonalBasin('zoneC', zone_C_region, 'Istrian Coast'),
        SimplePolygonalBasin('zoneD', zone_D_region, 'Depth Istria'),
        SimplePolygonalBasin('zoneE', zone_E_region, 'Quarnaro and North Dalmatia'),
        SimplePolygonalBasin('zoneF', zone_F_region, 'West Middle Adri'),
        SimplePolygonalBasin('zoneG', zone_G_region, 'Middle Middle Adri'),
        SimplePolygonalBasin('zoneH', zone_H_region, 'Zone H'),
        SimplePolygonalBasin('zoneI', zone_I_region, 'Zone I'),
        SimplePolygonalBasin('zoneJ', zone_J_region, 'Middle Dalmatia'),
        SimplePolygonalBasin('zoneK', zone_K_region, 'Zone K'),
        SimplePolygonalBasin('zoneL', zone_L_region, 'Zone L'),
        SimplePolygonalBasin('zoneM', zone_M_region, 'Zone M'),
    )

    return zones


def plot_zones():
    import matplotlib.pyplot as plt
    import cartopy
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    n_points = 2_000

    plt.figure(dpi=600, figsize=(16, 9))
    axes = plt.axes(projection=ccrs.PlateCarree())
    zones = generate_zones()
    for i, zone in enumerate(zones):
        zone.plot(
            lon_window=(OUTSIDE_EAST, OUTSIDE_WEST),
            lat_window=(DOMAIN_LIMIT_SOUTH - 0.05, OUTSIDE_NORTH),
            lon_points=n_points,
            lat_points=2 * n_points,
            color='C{}'.format(i),
            axes=axes,
            alpha=0.5
        )

    land = cfeature.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='10m',
        facecolor=(0.5, 0.5, 0.5, 0.8),
        edgecolor="black"
    )

    axes.add_feature(land, zorder=5)
    plt.savefig('zones.pdf')
    plt.close()


def get_zone_by_char(char_str):
    for zone in generate_zones():
        if zone.name.lower() == 'zone{}'.format(char_str.lower()):
            return zone
    raise ValueError('Zone {} not found'.format(char_str))


def generate_basins():
    basins = []

    def create_basin(bathymetric_region_name, zone):
        if bathymetric_region_name is None:
            basins.append(SimpleBasin(
                'z{}'.format(len(basins) + 1),
                zone.region,
                zone.extended_name
            ))
            return

        bathymetric_regions = generate_bathymetric_regions()
        bathymetric_region = bathymetric_regions[bathymetric_region_name]

        basin_name = '{} - {}'.format(zone.extended_name, bathymetric_region_name)

        basin_region = zone.region.intersect(bathymetric_region)
        basins.append(SimpleBasin(
            'z{}'.format(len(basins) + 1),
            basin_region,
            basin_name
        ))

    zoneB = get_zone_by_char('B')
    zoneB1 = get_zone_by_char('B1')
    zoneB2 = get_zone_by_char('B2')
    zoneB3 = get_zone_by_char('B3')
    zoneB4 = get_zone_by_char('B4')

    for depth in ('costal', 'shallow'):
        for zone in (zoneB4, zoneB3, zoneB2, zoneB1):
            create_basin(depth, zone)

    create_basin('mid', zoneB4)

    # Here we create basin number 10, which is inside three different zones
    basin_name = '{} - {}'.format(zoneB.extended_name, 'mid')
    bathymetric_region = generate_bathymetric_regions()['mid']

    basin_region = (zoneB1.region + zoneB2.region + zoneB3.region).intersect(
        bathymetric_region
    )
    basins.append(SimpleBasin(
        'z{}'.format(len(basins) + 1),
        basin_region,
        basin_name
    ))

    create_basin('deep', zoneB)

    zoneA = get_zone_by_char('A')
    for depth in ('costal', 'shallow'):
        create_basin(depth, zoneA)

    create_basin(None, get_zone_by_char('C'))
    create_basin(None, get_zone_by_char('D'))
    create_basin(None, get_zone_by_char('E'))

    for depth in ('costal', 'shallow', 'mid', 'deep'):
        create_basin(depth, get_zone_by_char('F'))

    create_basin(None, get_zone_by_char('G'))

    for depth in ('costal', 'shallow', 'mid', 'deep'):
        create_basin(depth, get_zone_by_char('H'))

    create_basin(None, get_zone_by_char('I'))
    create_basin(None, get_zone_by_char('J'))

    for depth in ('costal', 'shallow', 'mid', 'deep'):
        create_basin(depth, get_zone_by_char('K'))

    create_basin(None, get_zone_by_char('L'))
    create_basin(None, get_zone_by_char('M'))

    nad = ComposedBasin('nad', basins, 'North Adriatic (CADEAU)')

    return nad


def plot_basins():
    import matplotlib.pyplot as plt
    import cartopy
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    n_points = 2_000

    plt.figure(dpi=600, figsize=(16, 9))
    axes = plt.axes(projection=ccrs.PlateCarree())
    zones = generate_basins()
    for i, zone in enumerate(zones):
        zone.plot(
            lon_window=(OUTSIDE_EAST, OUTSIDE_WEST),
            lat_window=(DOMAIN_LIMIT_SOUTH - 0.05, OUTSIDE_NORTH),
            lon_points=n_points,
            lat_points=2 * n_points,
            color='C{}'.format(i),
            axes=axes,
            alpha=1
        )

    land = cfeature.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='10m',
        facecolor=(0.5, 0.5, 0.5, 0.8),
        edgecolor="black"
    )

    axes.add_feature(land, zorder=5)
    plt.savefig('basins.pdf')
    plt.close()


nad = generate_basins()


if __name__ == '__main__':
    # plot_bathymetric_regions()
    # plot_zones()
    # plot_basins()
    pass
