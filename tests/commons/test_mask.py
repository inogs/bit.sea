import numpy as np
import pytest

from bitsea.commons.mask import FILL_VALUE
from bitsea.commons.mask import Mask
from bitsea.commons.mask import MaskBathymetry
from bitsea.commons.mesh import Mesh


@pytest.fixture(scope="module", autouse=True)
def mesh():
    mesh_shape = (11, 19, 31)

    xlevels = np.arange(0, mesh_shape[2], dtype=np.float32)
    xlevels = (
        xlevels
        * np.linspace(1, 1.05, mesh_shape[1], dtype=np.float32)[:, np.newaxis]
    )

    y_levels = np.arange(0, mesh_shape[1], dtype=np.float32)
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], mesh_shape[1:])

    z_levels = np.geomspace(1, 1000, mesh_shape[0])

    mesh = Mesh.from_levels(xlevels, ylevels, z_levels)

    return mesh


@pytest.fixture(scope="module", autouse=True)
def mask(mesh):
    mask_array = np.zeros(mesh.shape, dtype=bool)
    for i in range(mesh.grid.shape[0]):
        for j in range(mesh.grid.shape[1]):
            if "3" in str(i) or "5" in str(i):
                if "7" in str(j) or "2" in str(j):
                    mask_array[:, i, j] = True

    for layer in range(1, mesh.shape[0]):
        mask_array[layer, :layer, :] = False
        mask_array[layer, -layer:, :] = False
        mask_array[layer, :, :layer] = False
        mask_array[layer, :, -layer:] = False

    return Mask(
        mesh.grid,
        zlevels=mesh.zlevels,
        mask_array=mask_array,
        allow_broadcast=True,
    )


def test_mask_follows_allow_broadcast(mesh):
    with pytest.raises(ValueError):
        Mask(
            mesh.grid,
            zlevels=mesh.zlevels,
            mask_array=[True],
            allow_broadcast=False,
        )

    mask_array = np.ones(shape=mesh.grid.shape, dtype=bool)
    mask = Mask(
        mesh.grid,
        zlevels=mesh.zlevels,
        mask_array=mask_array,
        allow_broadcast=True,
    )
    assert mask.ndim == 3


@pytest.mark.uses_test_data
def test_mask_from_file(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "nonregular_mask.nc"

    test_mask = Mask.from_file(mask_file)
    assert not test_mask.is_regular()


@pytest.mark.uses_test_data
def test_mask_from_file_when_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)
    assert test_mask.is_regular()
    assert isinstance(test_mask, Mask)


def test_mask_attributes(mask):
    assert len(mask.shape) == 3
    assert mask.shape == mask.as_array().shape
    assert mask.ndim == 3
    assert mask.dtype == bool


def test_mask_checks_zlevels(mesh):
    grid = mesh.grid
    zlevels = np.arange(10)[:, np.newaxis, np.newaxis]
    mesh_shape = (10,) + grid.shape
    zlevels_broadcast = np.broadcast_to(zlevels, mesh_shape)
    mask_array = np.random.choice([True, False], size=mesh_shape)

    Mask(grid, zlevels=zlevels[:, 0, 0], mask_array=mask_array)

    with pytest.raises(ValueError):
        Mask(grid, zlevels=zlevels_broadcast, mask_array=mask_array)


def convert_lon_lat_wetpoint_indices_reference_implementation(
    lon, lat, mask, maxradius=2
):
    """
    This was the previous implementation of the convert_lon_lat_wetpoint_indices
    method.
    We keep it as a reference, and we test that our current implementation is
    coherent with this one.
    """
    # Indexes of the input lon, lat
    lon = float(lon)
    lat = float(lat)

    jp, ip = mask.convert_lat_lon_to_indices(lon=lon, lat=lat)

    if mask[0, jp, ip]:
        return ip, jp

    # Matrixes of indexes of the Mask
    Ilist = np.arange(mask.shape[2])
    II = np.tile(Ilist, (mask.shape[1], 1))

    Jlist = np.arange(mask.shape[1]).T
    JJ = np.tile(Jlist, (mask.shape[2], 1)).T

    IImask = II[mask[0, :, :]]
    JJmask = JJ[mask[0, :, :]]

    # Find distances from wet points
    distind = np.sqrt((ip - IImask) ** 2 + (jp - JJmask) ** 2)
    if maxradius is None:
        maxradius = distind.min()
    indd = distind <= maxradius
    # Limit to distance < maxradius)
    ipnarr = IImask[indd]
    jpnarr = JJmask[indd]
    # Assign the nearest wet points
    if len(ipnarr) > 0:
        distmask = distind[indd]
        indmin = np.argmin(distmask)
        newip = ipnarr[indmin]
        newjp = jpnarr[indmin]
        return newip, newjp

    # If there aren't wet points with distance < maxradius, assign the non-wet point
    else:
        # print('WARNING: Using terrain point indexes, put maxradius=',
        #       distind.min(), ' or maxradius=None')
        return ip, jp


@pytest.mark.parametrize("max_radius", [0, 1, 2, 3, None])
@pytest.mark.filterwarnings("ignore:Around point")
def test_mask_convert_lon_lat_wetpoint_indices(mask, max_radius):
    lon_values = np.linspace(0, 30, 25)
    lat_values = np.linspace(0, 18, 23)[:, np.newaxis]

    lon_values, lat_values = np.broadcast_arrays(lon_values, lat_values)
    lon_values = lon_values.flatten()
    lat_values = lat_values.flatten()

    for lon, lat in zip(lon_values, lat_values):
        current_result = mask.convert_lon_lat_wetpoint_indices(
            lon=lon, lat=lat, max_radius=max_radius
        )
        previous_implementation = (
            convert_lon_lat_wetpoint_indices_reference_implementation(
                lon=lon, lat=lat, mask=mask, maxradius=max_radius
            )
        )

        assert current_result == previous_implementation


def test_mask_attribute_mask_is_deprecated(mask):
    with pytest.warns(DeprecationWarning):
        mask.mask


def test_mask_get_water_cells(mask):
    assert np.all(mask.get_water_cells() == mask)


def test_mask_get_sea_cells(mask):
    assert np.all(mask.get_sea_cells() == mask)


def test_mask_get_mask_at_level(mask):
    assert np.all(mask.mask_at_level(mask.zlevels[0] / 2) == mask[0])

    # Check that after the last level everything is False
    assert not np.any(mask.mask_at_level(mask.zlevels[-1] + 1))

    for z_i in np.linspace(3, 11, 10):
        i = mask.get_depth_index(z_i)
        assert np.all(mask.mask_at_level(z_i) == mask[i + 1, :])


def test_bathymetry_in_cells(mask):
    cell_bathymetry = mask.bathymetry_in_cells()
    for i in range(mask.shape[1]):
        for j in range(mask.shape[2]):
            if not np.any(mask[:, i, j]):
                assert cell_bathymetry[i, j] == 0
            else:
                assert np.count_nonzero(mask[:, i, j]) == cell_bathymetry[i, j]


def test_rough_bathymetry(mask):
    rough_bathymetry = mask.rough_bathymetry()
    for i in range(mask.shape[1]):
        for j in range(mask.shape[2]):
            if not np.any(mask[:, i, j]):
                current_depth = 0
            else:
                max_cell = np.where(mask[:, i, j])[0][-1]
                current_depth = mask.zlevels[max_cell]
            assert rough_bathymetry[i, j] == current_depth


def test_bathymetry(mask):
    bathymetry = mask.bathymetry()
    for i in range(mask.shape[1]):
        for j in range(mask.shape[2]):
            if not np.any(mask[:, i, j]):
                current_depth = FILL_VALUE
            else:
                current_depth = np.sum(mask.e3t[:, i, j], where=mask[:, i, j])
            assert bathymetry[i, j] == current_depth


def bathymetry_reference_implementation(mask):
    cells_bathy = mask.bathymetry_in_cells()
    _, jpj, jpi = mask.shape
    Bathy = np.ones((jpj, jpi), np.float32) * 1.0e20
    for ji in range(jpi):
        for jj in range(jpj):
            max_lev = cells_bathy[jj, ji]
            if max_lev > 0:
                Bathy[jj, ji] = mask.e3t[:max_lev, jj, ji].sum()
    return Bathy


def test_bathymetry_is_coherent_with_the_reference_implementation(mask):
    assert np.allclose(
        mask.bathymetry(), bathymetry_reference_implementation(mask)
    )


def test_cut_at_level(mask):
    for i in range(mask.shape[0]):
        layer = mask.cut_at_level(i)
        assert layer.depth == mask.zlevels[i]
        assert layer.thickness == mask.dz[i]
        assert np.all(layer[:] == mask[i])


def test_mask_coastline(mask):
    coastline = mask.coastline(depth=0)
    assert len(coastline) > 0


def test_bathymetry_from_mask(mask):
    bathymetry = MaskBathymetry(mask)

    for i in range(mask.grid.shape[0]):
        for j in range(mask.grid.shape[1]):
            lon, lat = mask.xlevels[i, j], mask.ylevels[i, j]

            assert bathymetry.is_inside_domain(lon=lon, lat=lat)

    assert not bathymetry.is_inside_domain(
        lon=np.max(mask.xlevels) + 1, lat=mask.ylevels[0, 0]
    )


def test_bathymetry_is_zero_on_land(mask):
    bathymetry = MaskBathymetry(mask)

    land_points_lon = mask.xlevels[~mask[0]]
    land_points_lat = mask.ylevels[~mask[0]]

    assert np.all(bathymetry(lon=land_points_lon, lat=land_points_lat) == 0)


def test_bathymetry_vectorizes(mask):
    bathymetry = MaskBathymetry(mask)

    min_lon = np.min(mask.xlevels)
    max_lon = np.max(mask.xlevels)

    min_lat = np.min(mask.ylevels)
    max_lat = np.max(mask.ylevels)

    lon_test = np.array(
        [min_lon, (max_lon - min_lon) / 3, (max_lon - min_lon) / 3 * 2, max_lon]
    )
    lat_test = np.array([min_lat, (max_lat - min_lat) / 2, max_lat])

    lon_broadcast, lat_broadcast = np.broadcast_arrays(
        lon_test, lat_test[:, np.newaxis]
    )

    bathymetry_values = bathymetry(lon=lon_broadcast, lat=lat_broadcast)
    for i in range(lat_test.shape[0]):
        for j in range(lon_test.shape[0]):
            assert (
                bathymetry(lon=lon_test[j], lat=lat_test[i])
                == bathymetry_values[i, j]
            )
