import datetime
import pickle

import netCDF4
import numpy as np

from bitsea.commons.Timelist import TimeList
from bitsea.validation.multirun.plot_profiles.tools.data_object import (
    NetCDFDataObject,
)
from bitsea.validation.multirun.plot_profiles.tools.data_object import (
    PickleDataObject,
)
from bitsea.validation.multirun.plot_profiles.tools.data_object import (
    get_data_object,
)


def test_get_data_object_routes_to_pickle(tmp_path):
    expected = np.ma.arange(2 * 2 * 3 * 4 * 5).reshape(2, 2, 3, 4, 5)
    time_steps = TimeList(
        [
            datetime.datetime(2024, 1, 1, 12, 0, 0),
            datetime.datetime(2024, 1, 2, 12, 0, 0),
        ],
        forceFrequency="daily",
    )

    with open(tmp_path / "N1p.pkl", "wb") as f:
        pickle.dump([expected, time_steps], f)

    data_object = get_data_object(tmp_path, "N1p")

    assert isinstance(data_object, PickleDataObject)
    np.testing.assert_array_equal(
        data_object.get_values(
            slice(None), slice(None), slice(None), slice(None), slice(None)
        ),
        expected,
    )
    assert data_object.get_time_steps() == time_steps.Timelist


def test_get_data_object_routes_to_netcdf(tmp_path):
    time_steps = [
        datetime.datetime(2024, 1, 1, 12, 0, 0),
        datetime.datetime(2024, 1, 2, 12, 0, 0),
    ]
    expected = []
    for index, time_step in enumerate(time_steps):
        values = np.arange(2 * 3 * 4 * 5, dtype=np.float32).reshape(
            2, 3, 4, 5
        ).copy()
        values += index * 1000
        expected.append(values)
        filename = tmp_path / time_step.strftime("ave.%Y%m%d-%H:%M:%S.stat_profiles.nc")
        with netCDF4.Dataset(filename, "w") as dataset:
            dataset.createDimension("sub", values.shape[0])
            dataset.createDimension("coast", values.shape[1])
            dataset.createDimension("z", values.shape[2])
            dataset.createDimension("stat", values.shape[3])

            variable = dataset.createVariable(
                "N1p", "f4", ("sub", "coast", "z", "stat")
            )
            variable[:] = np.ascontiguousarray(values)

    data_object = get_data_object(tmp_path, "N1p")

    assert isinstance(data_object, NetCDFDataObject)
    np.testing.assert_array_equal(
        data_object.get_values(
            slice(None), slice(None), slice(None), slice(None), slice(None)
        ),
        np.ma.stack(expected, axis=0),
    )
    assert data_object.get_time_steps() == time_steps
