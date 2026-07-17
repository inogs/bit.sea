from abc import ABC, abstractmethod
from collections.abc import Iterable
from contextlib import closing
import netCDF4
import numpy as np
from pathlib import Path
from typing import Union

from bitsea.commons.Timelist import TimeList
from bitsea.timeseries.plot import read_pickle_file


SLICE_TYPE = Union[slice, int]


class InvalidAxisSpecified(Exception):
    pass


class DataObject(ABC):

    @abstractmethod
    def get_values(self, time_steps: SLICE_TYPE, basin: SLICE_TYPE,
                   level_index: SLICE_TYPE, indicator: SLICE_TYPE,
                   coasts: SLICE_TYPE = 1):
        raise NotImplementedError

    @abstractmethod
    def get_time_steps(self):
        raise NotImplementedError

    @abstractmethod
    def is_2d(self) -> bool:
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return

    @staticmethod
    def get_axis(axis_name: str, while_fixing: Union[Iterable, None] = None):
        """
        Return the index of a specified axis identified by its name.

        :param axis_name: The name of the axis; one among "time", "basin",
        "coasts", "level" and "indicator"
        :param while_fixing: if the output has been previously sliced removing
        some indices, it is possible to insert here the name of the axis that
        have been sliced. In this way, the routine keeps track of the axes that
        have been removed and rescales the final output in order of obtaining
        the correct indices. So, for example, "level" has index 4; if "basin",
        which has index 2, has been fixed, then level is rescaled to index 3.

        :return: the index of the axis described by its name
        """
        axis_order = ('time', 'basin', 'coasts', 'level', 'indicator')

        if while_fixing is None:
            while_fixing = set()

        for current_axis_name in while_fixing:
            if current_axis_name not in axis_order:
                raise InvalidAxisSpecified(
                    'Invalid axis "{}" specified in the "with_fixed_axes" '
                    'argument. The only allowed values names are: {}'.format(
                        current_axis_name,
                        ', '.join(axis_order)
                    )
                )

        if axis_name not in axis_order:
            raise InvalidAxisSpecified(
                'Invalid axis "{}" specified. The only allowed values '
                'are: '.format(axis_name, ', '.join(axis_order))
            )

        axis_index = 0
        for name in axis_order:
            if name == axis_name:
                return axis_index
            if name not in while_fixing:
                axis_index += 1

        assert False, 'Internal error; there is a bug in the code!'


class PickleDataObject(DataObject):
    BFMv5_dict = {
        'Ac': 'ALK',
        'ppn': 'netPPYc',
        'ppg': 'ruPPYc',
        'ppb': 'ruPBAc',
        'CaCO3flux_dic': 'rcalCARc'
    }

    BFMv2_dict = {
        'ALK': 'Ac',
        'netPPYc': 'ppn',
        'netPPYc2': 'ppn',
        'ruPPYc': 'ppg',
        'ruPBAc': 'ppb',
        'rcalCARc': 'CaCO3flux_dic'
    }

    def __init__(self, dir_path, var_name):
        self._dir_path = Path(dir_path)
        self.var_name = var_name

        self._loaded = False
        self._data = None
        self._time_steps = None
        self._is_2d = None

    @classmethod
    def _iter_var_name_candidates(cls, var_name):
        yield var_name

        for alternative_names in (cls.BFMv2_dict, cls.BFMv5_dict):
            if var_name in alternative_names:
                yield alternative_names[var_name]

    @property
    def dir_path(self):
        return self._dir_path

    @dir_path.setter
    def dir_path(self, dir_path):
        self._dir_path = Path(dir_path)

    @classmethod
    def get_matching_filename(cls, dir_path, var_name):
        dir_path = Path(dir_path)
        for current_var_name in cls._iter_var_name_candidates(var_name):
            for extension in (".pkl", ".pickle"):
                candidate = dir_path / f"{current_var_name}{extension}"
                if candidate.exists():
                    return candidate
        return None

    def load(self):
        if self._loaded:
            return

        final_filename = self.get_matching_filename(self.dir_path, self.var_name)
        if final_filename is None:
            raise IOError(
                'No pickle file found for variable "{}" in directory "{}"'.format(
                    self.var_name, self.dir_path
                )
            )

        self._data, time_steps = read_pickle_file(final_filename)
        if not hasattr(self._data, 'mask'):
            self._data = np.ma.masked_invalid(self._data)
        self._time_steps = time_steps.Timelist
        self._loaded = True

    def get_values(self, time_steps, basin, level_index, indicator, coasts=1):
        if not self._loaded:
            self.load()
        return self._data[time_steps, basin, coasts, level_index, indicator]

    def get_time_steps(self):
        if not self._loaded:
            self.load()
        return self._time_steps

    def is_2d(self):
        if self._is_2d is None:
            if not self._loaded:
                self.load()
            self._is_2d = np.all(np.ma.getmaskarray(self._data)[0, :, :, 1:])
        return self._is_2d


class NetCDFDataObject(PickleDataObject):
    _search_patterns = ("*.stat_profiles.nc", "*.stat_profiles.netcdf")

    @classmethod
    def get_matching_filenames(cls, dir_path):
        dir_path = Path(dir_path)
        filenames = []
        for search_pattern in cls._search_patterns:
            filenames.extend(sorted(dir_path.glob(search_pattern)))
        return tuple(filenames)

    @classmethod
    def _read_time_steps(cls, dir_path):
        dir_path = str(Path(dir_path))
        filenames = []
        time_steps = []
        for search_pattern in cls._search_patterns:
            try:
                time_list = TimeList.fromfilenames(None, dir_path, search_pattern)
            except AssertionError:
                continue
            filenames.extend(time_list.filelist)
            time_steps.extend(time_list.Timelist)

        if not filenames:
            raise IOError(
                'No netCDF stat profile files found in directory "{}"'.format(
                    dir_path
                )
            )

        sorted_time_files = sorted(zip(time_steps, filenames), key=lambda x: x[0])
        time_steps, filenames = zip(*sorted_time_files)
        return TimeList(list(time_steps)), tuple(filenames)

    def _resolve_variable_name(self, filename):
        with closing(netCDF4.Dataset(filename, "r")) as dataset:
            for current_var_name in self._iter_var_name_candidates(self.var_name):
                if current_var_name in dataset.variables:
                    return current_var_name

        raise IOError(
            'Variable "{}" not found in netCDF stat profile file "{}"'.format(
                self.var_name, filename
            )
        )

    def load(self):
        if self._loaded:
            return

        time_steps, filenames = self._read_time_steps(self.dir_path)
        resolved_var_name = self._resolve_variable_name(filenames[0])

        data_list = []
        for filename in filenames:
            with closing(netCDF4.Dataset(filename, "r")) as dataset:
                if resolved_var_name not in dataset.variables:
                    raise IOError(
                        'Variable "{}" not found in netCDF stat profile file '
                        '"{}"'.format(resolved_var_name, filename)
                    )
                data_list.append(
                    np.ma.array(dataset.variables[resolved_var_name][:], copy=True)
                )

        self._data = np.ma.stack(data_list, axis=0)
        if not hasattr(self._data, 'mask'):
            self._data = np.ma.masked_invalid(self._data)
        self._time_steps = time_steps.Timelist
        self._loaded = True


def get_data_object(dir_path, var_name):
    if PickleDataObject.get_matching_filename(dir_path, var_name) is not None:
        return PickleDataObject(dir_path, var_name)

    if NetCDFDataObject.get_matching_filenames(dir_path):
        return NetCDFDataObject(dir_path, var_name)

    raise IOError(
        'No supported input file found for variable "{}" in directory "{}"'.format(
            var_name, dir_path
        )
    )
