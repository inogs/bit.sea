from abc import ABC, abstractmethod
from collections.abc import Iterable
import numpy as np
from pathlib import Path
from typing import Union

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

    @property
    def dir_path(self):
        return self._dir_path

    @dir_path.setter
    def dir_path(self, dir_path):
        self._dir_path = Path(dir_path)

    def load(self):
        if self._loaded:
            return

        filename_candidates = []
        main_filename = self.dir_path / (self.var_name + ".pkl")
        filename_candidates.append(main_filename)

        data_class = self.__class__
        for alternative_names in (data_class.BFMv2_dict, data_class.BFMv5_dict):
            if self.var_name in alternative_names:
                new_var_name = alternative_names[self.var_name]
                alternative_filename = self.dir_path / (new_var_name + ".pkl")
                filename_candidates.append(alternative_filename)
        for candidate in filename_candidates:
            if candidate.exists():
                final_filename = candidate
                break
        else:
            raise IOError('File "{}" not found'.format(main_filename))

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
