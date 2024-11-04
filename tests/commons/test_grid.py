import numpy as np
import pytest
from netCDF4 import Dataset

from bitsea.commons.grid import IrregularGrid
from bitsea.commons.grid import IrregularMaskLayer
from bitsea.commons.grid import MaskLayer
from bitsea.commons.grid import OutsideDomain
from bitsea.commons.grid import RegularGrid
from bitsea.commons.grid import RegularMaskLayer


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
        IrregularGrid(
            xlevels=xlevels, ylevels=np.asarray(ylevels, dtype=np.float64)
        )

    # Do not accept 1D arrays even if they broadcast (because they have the
    # shape
    with pytest.raises(ValueError):
        IrregularGrid(
            xlevels=np.linspace(0, 1, 10), ylevels=np.linspace(0, 1, 10)
        )


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
                lon=current_lon, lat=current_lat
            )

            assert p_coords == (i, j)


def test_convert_lon_lat_to_indices_array(grid):
    i_indices = np.array([0, 10, 42, 74])[:, np.newaxis]
    j_indices = np.array([[0, 10, 25, 99, 129]])

    lon_test = grid.xlevels[i_indices, j_indices] + 0.0002
    lat_test = grid.ylevels[i_indices, j_indices] - 0.0002

    p_coords = grid.convert_lon_lat_to_indices(lon=lon_test, lat=lat_test)

    assert np.all(p_coords[0] == i_indices)
    assert np.all(p_coords[1] == j_indices)


def test_convert_lon_lat_to_indices_regular(regular_grid):
    lat = np.linspace(0, 5, 100)

    lon = np.linspace(0, 6, 300)
    lon *= lon

    grid = RegularGrid(lon=lon, lat=lat)

    for i in (0, 10, 25, 99):
        current_lat = float(lat[i]) - 0.0002
        for j in (0, 10, 100, 299):
            current_lon = float(lon[j]) + 0.0002

            p_coords = grid.convert_lon_lat_to_indices(
                lon=current_lon, lat=current_lat
            )

            assert p_coords == (i, j)


def test_convert_lon_lat_to_indices_array_regular(regular_grid):
    i_values = np.array((0, 10, 25, 99), dtype=int)
    j_values = np.array((0, 10, 100, 299), dtype=int)
    lat_values = regular_grid.lat[i_values]
    lon_values = regular_grid.lon[j_values]

    p_coords = regular_grid.convert_lon_lat_to_indices(
        lon=lon_values, lat=lat_values
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
    lat = 4.0
    with pytest.raises(OutsideDomain):
        regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat)

    lon = 20.0
    lat = np.array((4, 5, 6, 7, 8), dtype=np.float32)
    with pytest.raises(OutsideDomain):
        regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat)


@pytest.mark.uses_test_data
def test_e1t_is_read_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    grid = IrregularGrid.from_file(mask_file)

    with Dataset(mask_file, "r") as f:
        e1t = f.variables["e1t"][0, 0].data

    assert np.allclose(e1t, grid.e1t, rtol=1e-5)


@pytest.mark.uses_test_data
def test_e2t_is_read_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    grid = IrregularGrid.from_file(mask_file)

    with Dataset(mask_file, "r") as f:
        e1t = f.variables["e2t"][0, 0].data

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
        lon=file_grid.lon, lat=file_grid.lat, e1t=None, e2t=None
    )
    assert np.allclose(file_grid.e1t, new_grid.e1t, rtol=1e-4)


@pytest.mark.uses_test_data
def test_e2t_can_be_computed_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    file_grid = RegularGrid.from_file(mask_file)
    new_grid = RegularGrid(
        lon=file_grid.lon, lat=file_grid.lat, e1t=None, e2t=None
    )
    assert np.allclose(file_grid.e2t, new_grid.e2t, rtol=1e-4)


def test_mask_layer_irregular(grid):
    mask = np.ones(grid.shape, dtype=bool)
    depth = 100.0
    thickness = 1.0
    mask_layer = MaskLayer.from_grid(
        grid, depth=depth, thickness=thickness, mask=mask
    )

    assert isinstance(mask_layer, IrregularMaskLayer)

    assert mask_layer.shape == mask_layer.as_array().shape
    assert mask_layer.depth == depth
    assert mask_layer.thickness == thickness


def test_mask_layer_regular(regular_grid):
    mask = np.ones(regular_grid.shape, dtype=bool)
    depth = 100.0
    thickness = 1.0
    mask_layer = MaskLayer.from_grid(
        regular_grid, depth=depth, thickness=thickness, mask=mask
    )

    assert isinstance(mask_layer, RegularMaskLayer)

    assert mask_layer.shape == mask_layer.as_array().shape
    assert mask_layer.depth == depth
    assert mask_layer.thickness == thickness


def test_mask_layer_checks_mask_shape(grid, regular_grid):
    mask = np.ones(tuple(i + 1 for i in grid.shape), dtype=bool)
    depth = 100.0
    thickness = 1.0
    with pytest.raises(ValueError):
        MaskLayer.from_grid(grid, depth=depth, thickness=thickness, mask=mask)

    regular_mask = np.ones(tuple(i + 1 for i in regular_grid.shape), dtype=bool)
    with pytest.raises(ValueError):
        MaskLayer.from_grid(
            regular_grid, depth=depth, thickness=thickness, mask=regular_mask
        )


def test_regular_mask_layer_checks_grid_is_regular(grid):
    mask = np.ones(grid.shape, dtype=bool)
    with pytest.raises(ValueError):
        RegularMaskLayer.from_grid(grid, depth=100.0, thickness=1.0, mask=mask)


def test_irregular_mask_layer_checks_grid_is_regular(regular_grid):
    mask = np.ones(regular_grid.shape, dtype=bool)
    with pytest.raises(ValueError):
        IrregularMaskLayer.from_grid(
            regular_grid, depth=100.0, thickness=1.0, mask=mask
        )


def test_regular_convert_lat_lon(regular_grid):
    lon_test = np.linspace(regular_grid.lon[0], regular_grid.lon[-1], 50)
    lat_test = np.linspace(regular_grid.lat[0], regular_grid.lat[-1], 100)

    for lon in lon_test:
        i, j = regular_grid.convert_lon_lat_to_indices(lon=lon, lat=lat_test[0])

        expected_j = np.argmin(np.abs(lon - regular_grid.lon))
        assert j == expected_j

    for lat in lat_test:
        i, j = regular_grid.convert_lon_lat_to_indices(lon=lon_test[0], lat=lat)

        expected_i = np.argmin(np.abs(lat - regular_grid.lat))
        assert i == expected_i


def test_regular_and_irregular_convert_i_j_to_lon_lat(regular_grid):
    irregular = IrregularGrid(
        xlevels=regular_grid.xlevels, ylevels=regular_grid.ylevels
    )

    i = np.arange(0, regular_grid.shape[0])[:, np.newaxis]
    j = np.arange(0, regular_grid.shape[1])[np.newaxis, :]

    assert np.all(
        irregular.convert_i_j_to_lon_lat(i, j)[0] == irregular.xlevels[i, j]
    )
    assert np.all(
        irregular.convert_i_j_to_lon_lat(i, j)[1] == irregular.ylevels[i, j]
    )
    assert np.all(
        regular_grid.convert_i_j_to_lon_lat(i, j)[0]
        == regular_grid.xlevels[i, j]
    )
    assert np.all(
        regular_grid.convert_i_j_to_lon_lat(i, j)[1]
        == regular_grid.ylevels[i, j]
    )


def test_regular_and_irregular_convert_lat_lon_to_i_j(regular_grid):
    irregular = IrregularGrid(
        xlevels=regular_grid.xlevels, ylevels=regular_grid.ylevels
    )

    lon_test = np.linspace(regular_grid.lon[0], regular_grid.lon[-1], 50)[
        np.newaxis, :
    ]
    lat_test = np.linspace(regular_grid.lat[0], regular_grid.lat[-1], 100)[
        :, np.newaxis
    ]

    t1 = regular_grid.convert_lon_lat_to_indices(lon=lon_test, lat=lat_test)
    t2 = irregular.convert_lon_lat_to_indices(lon=lon_test, lat=lat_test)

    for i in range(2):
        assert np.all(t1[i] == t2[i])
