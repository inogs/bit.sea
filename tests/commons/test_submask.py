import numpy as np
import pytest

from bitsea.basins.V2 import adr
from bitsea.basins.V2 import med
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask


@pytest.mark.uses_test_data
def test_mask_from_file_when_regular(test_data_dir):
    mask_dir = test_data_dir / "masks"
    mask_file = mask_dir / "regular_mask.nc"

    test_mask = Mask.from_file(mask_file)

    test_submask = SubMask(med, test_mask)
    assert np.all(test_mask[test_submask[:]])

    test_submask2 = SubMask(adr, test_mask)

    assert np.all(test_submask[test_submask2[:]])
    assert not np.all(test_submask2[test_submask[:]])
