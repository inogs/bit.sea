from pathlib import Path

import numpy as np
import pytest

from bitsea.basins.region import Polygon


@pytest.fixture
def quadrilateral():
    lon = [0, 1, 1, 0]
    lat = [-1, -1, 1, 0]

    return Polygon(lon, lat)


def test_polygon_single_value(quadrilateral):
    assert quadrilateral.is_inside(lon=0.5, lat=-0.1)
    assert not quadrilateral.is_inside(lon=0.5, lat=0.9)
    assert not quadrilateral.is_inside(lon=2, lat=3)


def test_polygon_vector(quadrilateral):
    lon_points = [0, 0.5, 0.5, -2, 0.9, 1]
    lat_points = [0.1, 0, -1.2, 0, 0.9, 2]

    test_values = quadrilateral.is_inside(lon=lon_points, lat=lat_points)

    assert np.all(test_values == [False, True, False, False, True, False])


def test_polygon_all_outside_rectangle(quadrilateral):
    lon_points = list(range(3, 20))
    lat = 0

    test_values = quadrilateral.is_inside(lon=lon_points, lat=lat)

    assert not np.any(test_values)


def test_polygon_matrix(quadrilateral):
    lon = np.linspace(-2, 2, 40)
    lat = np.linspace(-2, 2, 50)[:, np.newaxis]

    test_values = quadrilateral.is_inside(lon=lon, lat=lat)

    inside_square = np.logical_and(
        np.logical_and(lon >= 0, lon <= 1), np.logical_and(lat >= -1, lat <= 0)
    )

    inside_triangle = np.logical_and(
        np.logical_and(lon >= 0, lon <= 1),
        np.logical_and(lon - lat >= 0, lat >= 0),
    )
    expected_values = np.logical_or(inside_square, inside_triangle)

    assert np.all(test_values == expected_values)


def test_polygon_boundary():
    lon = [0, 3, 3, 2.3, 1, 0, 0]
    lat = [0, 0, 4, 7, 4.5, 6, 0]

    test_poly = Polygon(lon, lat)

    assert test_poly.border_longitudes == lon
    assert test_poly.border_latitudes == lat

    assert [p[0] for p in test_poly.borders] == lon
    assert [p[1] for p in test_poly.borders] == lat


@pytest.mark.uses_test_data
def test_read_wkt_file(test_data_dir: Path):
    wkt_file_path = test_data_dir / "wkt_polys.csv"
    with open(wkt_file_path, "r") as f:
        poly_dict = Polygon.read_WKT_file(f)
    assert len(poly_dict) == 2
    assert "Mambo" in poly_dict
    assert len(poly_dict["Mambo"].borders) == 8
