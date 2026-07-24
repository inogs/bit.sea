import json

import netCDF4
import numpy as np
import pytest

from bitsea.commons.grid import RegularGrid
from bitsea.commons.mask import MaskWithRivers


@pytest.fixture
def grid():
    lat = np.linspace(30, 45, 10)
    lon = np.linspace(-5, 30, 10)
    return RegularGrid(lat=lat, lon=lon)


@pytest.fixture
def zlevels():
    return np.array([1, 10, 100], dtype=float)


@pytest.fixture
def mask_data(grid, zlevels):
    shape = (len(zlevels),) + grid.shape
    mask = np.zeros(shape, dtype=bool)
    mask[:, 5, 5] = True  # Sea cell
    mask[:, 2, 2] = True  # River 1 cell
    mask[:, 3, 3] = True  # River 2 cell
    return mask


@pytest.fixture
def river_positions(grid):
    rivers = np.zeros(grid.shape, dtype=int)
    rivers[2, 2] = 1
    rivers[3, 3] = 2
    return rivers


def test_mask_with_rivers_initialization(
    grid, zlevels, mask_data, river_positions
):
    # Test without sources and names
    mask = MaskWithRivers(
        grid=grid,
        zlevels=zlevels,
        mask_array=mask_data,
        river_positions=river_positions,
    )
    assert not mask.has_river_sources()
    assert mask.n_rivers() == 2
    assert mask._river_names == {}


def test_mask_with_rivers_sources(grid, zlevels, mask_data, river_positions):
    river_sources = np.zeros(grid.shape, dtype=int)
    river_sources[2, 2] = 1

    mask = MaskWithRivers(
        grid=grid,
        zlevels=zlevels,
        mask_array=mask_data,
        river_positions=river_positions,
        river_sources=river_sources,
    )

    assert mask.has_river_sources()
    assert 1 in mask._sources
    assert mask._sources[1] == [(2, 2)]

    sources_map = mask.get_river_sources_map()
    assert np.all(sources_map == river_sources)


def test_mask_with_rivers_invalid_sources(
    grid, zlevels, mask_data, river_positions
):
    # Source not in river positions
    river_sources = np.zeros(grid.shape, dtype=int)
    river_sources[0, 0] = 1
    with pytest.raises(ValueError, match="not part of any river"):
        MaskWithRivers(
            grid=grid,
            zlevels=zlevels,
            mask_array=mask_data,
            river_positions=river_positions,
            river_sources=river_sources,
        )

    # Source ID mismatch
    river_sources = np.zeros(grid.shape, dtype=int)
    river_sources[2, 2] = 2  # Should be 1 according to river_positions[2,2]
    with pytest.raises(
        ValueError,
        match="Trying to set the source of river 2 on the stem of river 1",
    ):
        MaskWithRivers(
            grid=grid,
            zlevels=zlevels,
            mask_array=mask_data,
            river_positions=river_positions,
            river_sources=river_sources,
        )


def test_mask_with_rivers_names(grid, zlevels, mask_data, river_positions):
    river_names = {1: "Po", 2: "Rhone"}
    mask = MaskWithRivers(
        grid=grid,
        zlevels=zlevels,
        mask_array=mask_data,
        river_positions=river_positions,
        river_names=river_names,
    )
    assert mask._river_names == {1: "Po", 2: "Rhone"}


def test_mask_with_rivers_copy(grid, zlevels, mask_data, river_positions):
    river_sources = np.zeros(grid.shape, dtype=int)
    river_sources[2, 2] = 1
    river_names = {1: "Po", 2: "Rhone"}

    mask = MaskWithRivers(
        grid=grid,
        zlevels=zlevels,
        mask_array=mask_data,
        river_positions=river_positions,
        river_sources=river_sources,
        river_names=river_names,
    )

    copied = mask.copy()
    assert copied.has_river_sources()
    assert copied._sources == mask._sources
    assert copied._river_names == mask._river_names
    assert np.all(copied.get_river_sources_map() == river_sources)


def test_mask_with_rivers_netcdf_roundtrip(
    tmp_path, grid, zlevels, mask_data, river_positions
):
    river_sources = np.zeros(grid.shape, dtype=int)
    river_sources[2, 2] = 1
    river_sources[3, 3] = 2
    river_names = {1: "Po", 2: "Rhone"}

    mask = MaskWithRivers(
        grid=grid,
        zlevels=zlevels,
        mask_array=mask_data,
        river_positions=river_positions,
        river_sources=river_sources,
        river_names=river_names,
    )

    nc_path = tmp_path / "test_rivers.nc"
    mask.save_as_netcdf(nc_path)

    # Verify contents with netCDF4 directly
    with netCDF4.Dataset(nc_path, "r") as ds:
        assert "river_sources" in ds.variables
        assert np.all(ds.variables["river_sources"][:] == river_sources)
        assert ds.getncattr("river_names") == json.dumps(
            {"1": "Po", "2": "Rhone"}
        )

    # Load back using from_file
    loaded = MaskWithRivers.from_file(nc_path)
    assert loaded.has_river_sources()
    assert loaded._river_names == {1: "Po", 2: "Rhone"}
    assert np.all(loaded.get_river_sources_map() == river_sources)
    # Check that it also has the same sea mask and river positions
    assert np.all(loaded.get_water_cells() == mask.get_water_cells())
    assert np.all(loaded.as_array() == mask.as_array())


def test_get_river_sources_map_no_sources(
    grid, zlevels, mask_data, river_positions
):
    mask = MaskWithRivers(
        grid=grid,
        zlevels=zlevels,
        mask_array=mask_data,
        river_positions=river_positions,
    )
    with pytest.raises(
        ValueError, match="no information about the river sources"
    ):
        mask.get_river_sources_map()


def test_mask_with_rivers_from_file_no_sources(
    tmp_path, grid, zlevels, mask_data, river_positions
):
    # Create a NetCDF file WITHOUT river_sources variable
    nc_path = tmp_path / "test_no_sources.nc"

    with netCDF4.Dataset(nc_path, "w") as ds:
        ds.createDimension("z", len(zlevels))
        ds.createDimension("y", grid.shape[0])
        ds.createDimension("x", grid.shape[1])
        ds.createDimension("t", 1)

        v = ds.createVariable("nav_lev", "f4", ("z",))
        v[:] = zlevels
        v = ds.createVariable("nav_lat", "f4", ("y", "x"))
        v[:] = grid.ylevels
        v = ds.createVariable("nav_lon", "f4", ("y", "x"))
        v[:] = grid.xlevels
        v = ds.createVariable("tmask", "u1", ("z", "y", "x"))
        v[:] = mask_data
        v = ds.createVariable("rivers", "i4", ("y", "x"))
        v[:] = river_positions
        # NO river_sources variable here

    loaded = MaskWithRivers.from_file(nc_path, read_e3t=False)
    assert not loaded.has_river_sources()
