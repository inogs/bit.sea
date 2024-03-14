from abc import ABC, abstractmethod
from typing import Sequence
from typing import Union

from plot_inputs import PlotInputData


class InvalidFilterDescription(Exception):
    pass


class FilteredObject(PlotInputData, ABC):
    def __init__(self, input_data: PlotInputData):
        super().__init__(input_data._data_object)
        self._input_data = input_data

    def __enter__(self):
        self._input_data.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._input_data.__exit__(exc_type, exc_val, exc_tb)
        return


class Filter(ABC):
    @abstractmethod
    def get_filtered_object(self, input_data) -> FilteredObject:
        raise NotImplementedError


class SimpleFilter(Filter):
    @staticmethod
    @abstractmethod
    def initialize_from_string(args_str: Union[str, None]):
        raise NotImplementedError


class ComposedFilter(Filter):
    def __init__(self, filters: Sequence[Filter]):
        self._filters = tuple(filters)

    def get_filtered_object(self, input_data):
        current_data = input_data
        for current_filter in reversed(self._filters):
            current_data = current_filter.get_filtered_object(current_data)
        return current_data
