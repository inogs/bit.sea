from netCDF4 import Dataset
import numpy as np
import pytest

from bitsea.commons.grid import IrregularGrid
from bitsea.commons.grid import GeoPySidesCalculator
from bitsea.commons.grid import NemoGridSidesCalculator
from bitsea.commons.grid import OutsideDomain
from bitsea.commons.grid import RegularGrid
from bitsea.commons.grid import extend_from_average


@pytest.fixture
def grid():
    grid_shape = (75, 130)

    xlevels = np.linspace(-5, 5, grid_shape[1])
    xlevels = xlevels * np.linspace(1, 2, grid_shape[0])[:, np.newaxis]

    y_levels = np.linspace(20, 30, grid_shape[0])
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], grid_shape)

    return IrregularGrid(xlevels=xlevels, ylevels=ylevels)


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

    grid = IrregularGrid(xlevels=xlevels, ylevels=ylevels)

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

def test_init_with_wrongs_arguments():
    grid_shape = (19, 31)
    xlevels = np.linspace(0, 1, grid_shape[1], dtype=np.float32)
    xlevels = np.broadcast_to(xlevels, grid_shape)

    ylevels = np.linspace(0, 1, grid_shape[0], dtype=np.float32)
    ylevels = np.broadcast_to(ylevels[:, np.newaxis], grid_shape)

    # Must have the right shape
    with pytest.raises(ValueError):
        IrregularGrid(xlevels=xlevels[1:-1], ylevels=ylevels)

    # Must have same type
    with pytest.raises(ValueError):
        IrregularGrid(xlevels=xlevels, ylevels=np.asarray(ylevels, dtype=np.float64))

    # Do not accept 1D arrays even if they broadcast (because they have the
    # shape
    with pytest.raises(ValueError):
        IrregularGrid(xlevels=np.linspace(0, 1, 10), ylevels=np.linspace(0, 1, 10))


def test_init_grid_does_broadcast():
    grid_shape = (19, 31)
    xlevels = np.linspace(0, 1, grid_shape[1], dtype=np.float32)

    ylevels = np.linspace(0, 1, grid_shape[0], dtype=np.float32)
    ylevels = np.broadcast_to(ylevels[:, np.newaxis], grid_shape)

    IrregularGrid(xlevels=xlevels, ylevels=ylevels)


@pytest.mark.uses_test_data
def test_grid_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    grid = IrregularGrid.from_file(mask_file)
    assert not grid.is_regular()


@pytest.mark.uses_test_data
def test_regular_grid_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    grid = IrregularGrid.from_file(mask_file)
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


@pytest.mark.uses_test_data
def test_e1t_is_read_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    grid = IrregularGrid.from_file(mask_file)

    with Dataset(mask_file, 'r') as f:
        e1t = f.variables['e1t'][0, 0].data

    assert np.allclose(e1t, grid.e1t, rtol=1e-5)


@pytest.mark.uses_test_data
def test_e2t_is_read_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    grid = IrregularGrid.from_file(mask_file)

    with Dataset(mask_file, 'r') as f:
        e1t = f.variables['e2t'][0, 0].data

    assert np.allclose(e1t, grid.e2t, rtol=1e-5)


def test_e1t_can_be_computed(grid):
    assert grid.e1t is not None
    assert grid.e2t is not None


@pytest.mark.uses_test_data
def test_e1t_can_be_computed_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    file_grid = RegularGrid.from_file(mask_file)
    new_grid = RegularGrid(
        lon=file_grid.lon,
        lat=file_grid.lat,
        e1t=None,
        e2t=None
    )
    assert np.allclose(file_grid.e1t, new_grid.e1t, rtol=1e-4)


@pytest.mark.uses_test_data
def test_e2t_can_be_computed_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    file_grid = RegularGrid.from_file(mask_file)
    new_grid = RegularGrid(
        lon=file_grid.lon,
        lat=file_grid.lat,
        e1t=None,
        e2t=None
    )
    assert np.allclose(file_grid.e2t, new_grid.e2t, rtol=1e-4)


def test_different_algorithms_give_coherent_results(grid):
    c1 = NemoGridSidesCalculator()
    c2 = GeoPySidesCalculator(geodesic=False)
    c3 = GeoPySidesCalculator(geodesic=True)

    e1t_c1, e2t_c1 = c1(grid.xlevels, grid.ylevels)
    e1t_c2, e2t_c2  = c2(grid.xlevels, grid.ylevels)
    e1t_c3, e2t_c3 = c3(grid.xlevels, grid.ylevels)

    assert np.allclose(e1t_c1, e1t_c2, rtol=1e-2)
    assert np.allclose(e1t_c2, e1t_c3, rtol=1e-2)

    assert np.allclose(e2t_c1, e2t_c2, rtol=1e-1)
    assert np.allclose(e2t_c2, e2t_c3, rtol=1e-2)


def test_extend_from_average():
    t_array = np.linspace(0, 10, 7)
    t_array *= t_array

    v_array = extend_from_average(t_array)
    v_mean = (v_array[1:] + v_array[:-1]) / 2

    assert np.allclose(v_mean, t_array)
