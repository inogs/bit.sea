import os
from pathlib import Path

import numpy as np
import xarray as xr

from bitsea.commons.grid import RegularGrid
from bitsea.mhelpers.linear_shift import linear_shift

default_clim_n3n_nc = "/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/Nutrients/MedBGC-ins/Climatology_N3n_600_800m.nc"

clim_file = Path(os.getenv("CLIM_MED_BGC_INS_FILE", default_clim_n3n_nc))
if not clim_file.exists():
    print(clim_file)
    raise ValueError(
        "Environment variable CLIM_MED_BGC_INS_FILE must be defined"
    )


nc = xr.open_dataset(clim_file)
TheMask = RegularGrid(lon=nc.Lon, lat=nc.Lat)

N3n = nc.N3n.data


def nitrate_correction(p, Pres, Value):
    """
    Calculates the nitrate value of MedBGCins and corrects all the profile by the bottom value

    Bottom value is defined as mean in 600-800m
    Argument:
    * p * float profile object
    * Pres  * 1d ndarray
    * Value * 1d ndarray

    Return:
    * Pres  * numpy 1d array
    * Value * corrected values
    * Qc    * array of integers, 8
    * shift * float value
    """
    ji, jj = TheMask.convert_lon_lat_to_indices(p.lon, p.lat)

    N3n_600_800 = N3n[jj, ji]
    assert not np.isnan(N3n_600_800)

    ii = (Pres >= 600) & (Pres <= 800)
    if ii.sum() == 0:
        ii = Pres > 600
    N_Float_bottom = Value[ii].mean()

    shift = N_Float_bottom - N3n_600_800

    New_profile = linear_shift(Value, Pres, shift)
    Nqc = np.ones((len(Pres)), int) * 8
    return Pres, New_profile, Nqc, shift
