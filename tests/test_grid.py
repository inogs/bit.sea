import numpy as np

from bitsea.commons.mask import Grid
from bitsea.commons.mask import RegularGrid


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


def test_convert_lon_lat_to_indices():
    grid_shape = (75, 130)

    xlevels = np.linspace(-5, 5, grid_shape[1])
    xlevels = xlevels * np.linspace(1, 2, grid_shape[0])[:, np.newaxis]

    y_levels = np.linspace(20, 30, grid_shape[0])
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], grid_shape)

    grid = Grid(xlevels=xlevels, ylevels=ylevels)

    for i in (0, 10, 42, 74):
        for j in (0, 10, 25, 99, 129):
            current_lon = float(xlevels[i, j]) + 0.0002
            current_lat = float(ylevels[i, j]) - 0.0002

            p_coords = grid.convert_lon_lat_to_indices(
                lon=current_lon,
                lat=current_lat
            )

            assert p_coords == (i, j)


def test_convert_lon_lat_to_indices_array():
    grid_shape = (75, 130)

    xlevels = np.linspace(-5, 5, grid_shape[1])
    xlevels = xlevels * np.linspace(1, 2, grid_shape[0])[:, np.newaxis]

    y_levels = np.linspace(20, 30, grid_shape[0])
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], grid_shape)

    grid = Grid(xlevels=xlevels, ylevels=ylevels)

    i_indices = np.array([0, 10, 42, 74])[:, np.newaxis]
    j_indices = np.array([[0, 10, 25, 99, 129]])

    lon_test = xlevels[i_indices, j_indices] + 0.0002
    lat_test = ylevels[i_indices, j_indices] - 0.0002

    p_coords = grid.convert_lon_lat_to_indices(
        lon=lon_test,
        lat=lat_test
    )

    assert np.all(p_coords[0] == i_indices)
    assert np.all(p_coords[1] == j_indices)


def test_convert_lon_lat_to_indices_regular():
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


def test_convert_lon_lat_to_indices_array_regular():
    lat = np.linspace(0, 5, 100)

    lon = np.linspace(0, 6, 300)
    lon *= lon

    grid = RegularGrid(lon=lon, lat=lat)

    i_values = np.array((0, 10, 100, 299), dtype=int)
    j_values = np.array((0, 10, 25, 99), dtype=int)
    lon_values = lon[i_values]
    lat_values = lat[j_values]

    p_coords = grid.convert_lon_lat_to_indices(
        lon=lon_values,
        lat=lat_values
    )

    assert np.all(p_coords[0] == i_values)
    assert np.all(p_coords[1] == j_values)
