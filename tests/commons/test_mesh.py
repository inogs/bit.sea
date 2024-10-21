import netCDF4
import pytest
import numpy as np

from bitsea.commons.grid import IrregularGrid
from bitsea.commons.grid import Grid
from bitsea.commons.mesh import Mesh
from bitsea.commons.mesh import RegularMesh


@pytest.fixture
def grid():
    grid_shape = (75, 130)

    xlevels = np.linspace(-5, 5, grid_shape[1])
    xlevels = xlevels * np.linspace(1, 2, grid_shape[0])[:, np.newaxis]

    y_levels = np.linspace(20, 30, grid_shape[0])
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], grid_shape)

    return IrregularGrid(xlevels=xlevels, ylevels=ylevels)


@pytest.fixture
def mesh(grid):
    return Mesh(grid, 5. + np.arange(10) * 10.)


def test_mesh_grid_descriptor_interface(mesh):
    grid = mesh.grid

    for m in Grid.__abstractmethods__:
        m_method = getattr(Grid, m)
        if isinstance(m_method, property):
            if m == 'shape':
                assert grid.shape == mesh.shape[1:]
                continue
            if m == 'coordinate_dtype':
                assert grid.coordinate_dtype == mesh.coordinate_dtype
                continue
            assert np.allclose(getattr(grid, m), getattr(mesh, m))

    assert grid.is_regular() == mesh.is_regular()


def test_mesh_zlevels_must_be_non_negative(grid):
    zlevels = np.arange(10, 30, dtype=np.float32)
    zlevels[0] *= -1
    with pytest.raises(ValueError):
        Mesh(grid, zlevels)


def test_mesh_zlevels_must_be_strictly_increasing(grid):
    zlevels = np.arange(10, 30, dtype=np.float32)
    zlevels[5] = 13.8
    with pytest.raises(ValueError):
        Mesh(grid, zlevels)


def test_mesh_zlevels_must_be_1d(grid):
    zlevels = np.arange(10, 30, dtype=np.float32)
    zlevels = zlevels.reshape((-1, 1))
    with pytest.raises(ValueError):
        Mesh(grid, zlevels)


def test_mesh_computes_e3t_correctly(mesh):
    cell_heights = np.zeros(
        (mesh.zlevels.shape[0] + 1,), dtype=mesh.zlevels.dtype
    )
    cell_heights[1:] = np.cumsum(mesh.e3t[:, 0, 0])

    expected_zlevels = (cell_heights[1:] + cell_heights[:-1]) / 2.

    assert np.allclose(mesh.zlevels, expected_zlevels)


def test_mesh_can_handle_1d_e3t_arrays(grid):
    e3t = np.geomspace(1, 100, 9, dtype=np.float32)
    zlevels = e3t + 0.5

    mesh = Mesh(grid, zlevels, e3t=e3t)

    assert mesh.e3t.ndim == 3
    assert np.allclose(e3t[:, np.newaxis, np.newaxis], mesh.e3t)


def test_mesh_checks_if_e3t_has_the_right_shape(grid):
    e3t = np.ones((10, ) + grid.shape, dtype=np.float32)
    zlevels = np.arange(9) + 0.5

    with pytest.raises(ValueError):
        Mesh(grid, zlevels, e3t=e3t)

    e3t = np.ones((9, 11, 11), dtype=np.float32)
    with pytest.raises(ValueError):
        Mesh(grid, zlevels, e3t=e3t)


def test_mesh_dz_is_coherent(mesh):
    assert np.allclose(mesh.dz[:, np.newaxis, np.newaxis], mesh.e3t)


def test_mesh_get_depth_index(mesh):
    min_depth = mesh.zlevels[0]
    max_depth = mesh.zlevels[-1]

    for depth in np.linspace(min_depth, max_depth + 10., 50):
        depth_index = mesh.get_depth_index(depth)

        assert depth >= mesh.zlevels[depth_index]

        if depth_index != mesh.zlevels.shape[0] - 1:
            assert depth < mesh.zlevels[depth_index + 1]

    for depth in np.linspace(0, min_depth, 10):
        depth_index = mesh.get_depth_index(depth)

        assert depth_index == 0


def test_mesh_get_depth_index_vectorize(mesh):
    min_depth = mesh.zlevels[0]
    max_depth = mesh.zlevels[-1]

    d = np.linspace(min_depth, max_depth + 10., 50)

    d_indices = mesh.get_depth_index(d)

    for depth, depth_index in zip(d, d_indices):
        depth_index_standard = mesh.get_depth_index(depth)
        assert depth_index == depth_index_standard


def test_mesh_from_levels(mesh):
    new_mesh = Mesh.from_levels(mesh.xlevels, mesh.ylevels, mesh.zlevels)

    assert np.allclose(new_mesh.xlevels, mesh.xlevels)
    assert np.allclose(new_mesh.ylevels, mesh.ylevels)
    assert np.allclose(new_mesh.zlevels, mesh.zlevels)


@pytest.mark.uses_test_data
def test_mesh_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    file_mesh = Mesh.from_file(mask_file)
    file_grid = IrregularGrid.from_file(mask_file)

    assert np.allclose(file_mesh.xlevels, file_grid.xlevels)
    assert np.allclose(file_mesh.ylevels, file_grid.ylevels)

    assert not file_mesh.is_regular()


@pytest.mark.uses_test_data
def test_mesh_from_file_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    file_mesh = Mesh.from_file(mask_file)

    assert file_mesh.is_regular()


@pytest.mark.uses_test_data
def test_regular_mesh_checks_if_file_is_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    with pytest.raises(ValueError):
        RegularMesh.from_file(mask_file)


@pytest.mark.uses_test_data
def test_mesh_from_file_reading_e3t(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    file_mesh = Mesh.from_file(mask_file, read_e3t=True)

    with netCDF4.Dataset(mask_file, "r") as f:
        e3t = np.array(f.variables["e3t"][:], dtype=np.float32)

    assert np.allclose(e3t, file_mesh.e3t)


def test_regular_mesh_from_coordinates():
    lon = np.linspace(-5, 15, 50)
    lat = np.linspace(35, 46, 30)

    zlevels = np.geomspace(1, 200, 50)

    regular_mesh = RegularMesh.from_coordinates(
        lon=lon, lat=lat, zlevels=zlevels
    )

    assert np.allclose(regular_mesh.lon, lon)
    assert np.allclose(regular_mesh.lat, lat)
    assert np.allclose(regular_mesh.zlevels, zlevels)
