from pathlib import Path

import pytest

from bitsea.commons import netcdf4
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.interpolators import compose_methods
from bitsea.commons.interpolators import regular
from bitsea.commons.interpolators import space_interpolator_griddata
from bitsea.commons.mask import Mask


@pytest.mark.uses_test_data
def test_interpolator(tmp_path: Path, test_data_dir: Path):
    masks_dir = test_data_dir / "masks"
    Mask1 = Mask.from_file(masks_dir / "mask_006_014_reduced.nc")
    Mask2 = Mask.from_file(masks_dir / "interpolator_mask.nc")

    filename = test_data_dir / "ave.20241027-12:00:00.N1p.nc"

    VAR = DataExtractor(Mask1, filename, "N1p").values

    A = space_interpolator_griddata(Mask2, Mask1, VAR)
    B = regular(Mask1, Mask2, VAR, method="nearest")
    C = compose_methods(Mask1, Mask2, VAR)
    A[~Mask2] = 1.0e20

    netcdf4.write_3d_file(
        A, "N1p", tmp_path / "griddata.nc", Mask2, fillValue=1e20
    )
    netcdf4.write_3d_file(
        B, "N1p", tmp_path / "regular.nc", Mask2, fillValue=1e20
    )
    netcdf4.write_3d_file(
        C, "N1p", tmp_path / "composed.nc", Mask2, fillValue=1e20
    )
