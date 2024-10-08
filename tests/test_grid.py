import numpy as np
import pytest

from bitsea.commons.grid import Grid, OutsideDomain
from bitsea.commons.grid import RegularGrid


@pytest.fixture
def grid():
    grid_shape = (75, 130)

    xlevels = np.linspace(-5, 5, grid_shape[1])
    xlevels = xlevels * np.linspace(1, 2, grid_shape[0])[:, np.newaxis]

    y_levels = np.linspace(20, 30, grid_shape[0])
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], grid_shape)

    return Grid(xlevels=xlevels, ylevels=ylevels)


@pytest.fixture
def regular_grid():
    lat = np.linspace(0, 5, 100)
    lon = np.linspace(0, 6, 300)
    lon *= lon

    return RegularGrid(lon=lon, lat=lat)


def test_grid_init():
    grid_shape = (75, 130)

    xlevels = np.linspace(0, 1, grid_shape[1])
    xlevels = np.broadcast_to(xlevels, grid_shape)

    ylevels = np.linspace(0, 1, grid_shape[0])
    ylevels = np.broadcast_to(ylevels[:, np.newaxis], grid_shape)

    grid = Grid(xlevels=xlevels, ylevels=ylevels)

    assert np.allclose(xlevels, grid.xlevels)
    assert np.allclose(ylevels, grid.ylevels)

    assert grid.shape == grid_shape


def test_regular_grid_init():
    grid_shape = (75, 130)

    lat = np.linspace(0, 1, grid_shape[0])
    lon = np.linspace(0, 1, grid_shape[1])

    grid = RegularGrid(lon=lon, lat=lat)

    assert np.allclose(lon, grid.lon)
    assert np.allclose(lat, grid.lat)

    assert np.allclose(lon, grid.xlevels)
    assert np.allclose(lat[:, np.newaxis], grid.ylevels)

    assert grid.shape == grid_shape


def test_grid_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    grid = Grid.from_file(mask_file)
    assert not grid.is_regular()


def test_regular_grid_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    grid = Grid.from_file(mask_file)
    assert grid.is_regular()


def test_convert_lon_lat_to_indices(grid):
    for i in (0, 10, 42, 74):
        for j in (0, 10, 25, 99, 129):
            current_lon = float(grid.xlevels[i, j]) + 0.0002
            current_lat = float(grid.ylevels[i, j]) - 0.0002

            p_coords = grid.convert_lon_lat_to_indices(
                lon=current_lon,
                lat=current_lat
            )

            assert p_coords == (i, j)


def test_convert_lon_lat_to_indices_array(grid):
    i_indices = np.array([0, 10, 42, 74])[:, np.newaxis]
    j_indices = np.array([[0, 10, 25, 99, 129]])

    lon_test = grid.xlevels[i_indices, j_indices] + 0.0002
    lat_test = grid.ylevels[i_indices, j_indices] - 0.0002

    p_coords = grid.convert_lon_lat_to_indices(
        lon=lon_test,
        lat=lat_test
    )

    assert np.all(p_coords[0] == i_indices)
    assert np.all(p_coords[1] == j_indices)


def test_convert_lon_lat_to_indices_regular(regular_grid):
    lat = np.linspace(0, 5, 100)

    lon = np.linspace(0, 6, 300)
    lon *= lon

    grid = RegularGrid(lon=lon, lat=lat)

    for i in (0, 10, 100, 299):
        current_lon = float(lon[i]) + 0.0002
        for j in (0, 10, 25, 99):
            current_lat = float(lat[j]) - 0.0002

            p_coords = grid.convert_lon_lat_to_indices(
                lon=current_lon,
                lat=current_lat
            )

            assert p_coords == (i, j)


def test_convert_lon_lat_to_indices_array_regular(regular_grid):
    i_values = np.array((0, 10, 100, 299), dtype=int)
    j_values = np.array((0, 10, 25, 99), dtype=int)
    lon_values = regular_grid.lon[i_values]
    lat_values = regular_grid.lat[j_values]

    p_coords = regular_grid.convert_lon_lat_to_indices(
        lon=lon_values,
        lat=lat_values
    )

    assert np.all(p_coords[0] == i_values)
    assert np.all(p_coords[1] == j_values)


def test_outside_domain_lon_lat(grid):
    lon = -11
    lat = 21
    with pytest.raises(OutsideDomain):
        grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    lon = -5
    lat = 19
    with pytest.raises(OutsideDomain):
        grid.convert_lon_lat_to_indices(lon=lon, lat=lat)


def test_outside_domain_lon_lat_regular(regular_grid):
    lon = 38
    lat = 4
    with pytest.raises(OutsideDomain):
        regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    lon = 20
    lat = 8
    with pytest.raises(OutsideDomain):
        regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat)


def test_outside_domain_lon_lat_array(grid):
    lon = np.array((1, 2, 3, -11, 10), dtype=np.float32)
    lat = np.array((21, 21, 21, 21, 21), dtype=np.float32)
    with pytest.raises(OutsideDomain):
        grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    lon = np.array((-5, -5, -5, -5, -5), dtype=np.float32)
    lat = np.array((1, 10, 20, 50, 100), dtype=np.float32)
    with pytest.raises(OutsideDomain):
        grid.convert_lon_lat_to_indices(lon=lon, lat=lat)


def test_outside_domain_lon_lat_array_regular(regular_grid):
    lon = np.array((35, 36, 37, 38, 30), dtype=np.float32)
    lat = 4.
    with pytest.raises(OutsideDomain):
        regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    lon = 20.
    lat = np.array((4, 5, 6, 7, 8), dtype=np.float32)
    with pytest.raises(OutsideDomain):
        regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat)