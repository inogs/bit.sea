import numpy as np
import pytest

from bitsea.basins.region import Rectangle
from bitsea.basins.V2 import adr
from bitsea.basins.V2 import med
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask


@pytest.mark.uses_test_data
def test_submask_nesting(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    test_submask = SubMask(med, test_mask)
    assert np.all(test_mask[test_submask[:]])

    test_submask2 = SubMask(adr, test_mask)

    assert np.all(test_submask[test_submask2[:]])
    assert not np.all(test_submask2[test_submask[:]])


@pytest.mark.uses_test_data
def test_submask_sea_and_water_cells(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    test_submask = SubMask(adr, test_mask)

    assert np.all(test_submask.get_water_cells() == test_mask[:])
    assert np.all(test_submask.get_sea_cells() == test_mask[:])


@pytest.mark.uses_test_data
def test_submask_basin(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    test_submask = SubMask(adr, test_mask)

    assert test_submask.basin.get_uuid() == adr.get_uuid()


@pytest.mark.uses_test_data
def test_submask_from_square_cutting(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    mask_min_lon = np.min(test_mask.xlevels)
    mask_max_lon = np.max(test_mask.xlevels)

    mask_min_lat = np.min(test_mask.ylevels)
    mask_max_lat = np.max(test_mask.ylevels)

    for degrees in (1, 0.3, 0.1):
        n_x_squares = int((mask_max_lon - mask_min_lon) / degrees) + 1
        n_y_squares = int((mask_max_lat - mask_min_lat) / degrees) + 1

        rectangles = SubMask.from_square_cutting(test_mask, degrees)
        assert len(rectangles) == n_x_squares * n_y_squares

        assert np.abs(rectangles[0].latmin - mask_min_lat) < 1e-5
        assert np.abs(rectangles[0].lonmin - mask_min_lon) < 1e-5

        for rectangle in rectangles:
            assert np.abs(rectangle.latmax - rectangle.latmin - degrees) < 1e-5
            assert np.abs(rectangle.lonmax - rectangle.lonmin - degrees) < 1e-5


@pytest.mark.uses_test_data
def test_submask_from_square_cutting_checks_degrees(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    with pytest.raises(ValueError):
        SubMask.from_square_cutting(test_mask, -1)

    with pytest.raises(ValueError):
        SubMask.from_square_cutting(test_mask, 100)


@pytest.mark.uses_test_data
def test_submask_build_partition(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    mask_min_lon = np.min(test_mask.xlevels)
    mask_max_lon = np.max(test_mask.xlevels)
    lon_size = mask_max_lon - mask_min_lon

    mask_min_lat = np.min(test_mask.ylevels)
    mask_max_lat = np.max(test_mask.ylevels)

    r1 = Rectangle(
        lonmin=mask_min_lon,
        lonmax=mask_min_lon + 2 * lon_size / 3,
        latmin=mask_min_lat,
        latmax=mask_max_lat,
    )
    r2 = Rectangle(
        lonmin=mask_min_lon + lon_size / 3,
        lonmax=mask_max_lon,
        latmin=mask_min_lat,
        latmax=mask_max_lat,
    )

    sub_mask1 = SubMask(r1, test_mask)
    sub_mask2 = SubMask(r2, test_mask)

    partition = SubMask.build_partition(test_mask, [r1, r2])

    assert np.all(sub_mask1[:] == partition[0])
    assert np.all(sub_mask2[partition[1]])

    assert np.any(np.logical_and(sub_mask1[:], sub_mask2[:]))
    assert not np.any(np.logical_and(partition[0][:], partition[1][:]))
