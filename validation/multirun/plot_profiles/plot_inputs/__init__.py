from abc import ABC, abstractmethod
from typing import Union

from ..tools.data_object import DataObject


SLICE_TYPE = Union[slice, int]


class PlotInputData(ABC):
    def __init__(self, data_object: DataObject):
        self._data_object: DataObject = data_object

    @abstractmethod
    def get_time_steps(self):
        raise NotImplementedError

    @abstractmethod
    def get_values(self, time_steps: SLICE_TYPE, basin: SLICE_TYPE,
                   level_index: SLICE_TYPE):
        raise NotImplementedError

    def is_2d(self) -> bool:
        return self._data_object.is_2d()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return


