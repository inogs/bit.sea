from abc import ABC, abstractmethod
from typing import Sequence
from typing import Union

from tools.data_object import DataObject


class InvalidFilterDescription(Exception):
    pass


class FilteredObject(DataObject, ABC):
    def __init__(self, data_object: DataObject):
        self._original_data = data_object


class Filter(ABC):
    @abstractmethod
    def get_filtered_object(self, data_object) -> FilteredObject:
        raise NotImplementedError


class SimpleFilter(Filter):
    @staticmethod
    @abstractmethod
    def initialize_from_string(args_str: Union[str, None]):
        raise NotImplementedError


class ComposedFilter(Filter):
    def __init__(self, filters: Sequence[Filter]):
        self._filters = tuple(filters)

    def get_filtered_object(self, data_object):
        current_data = data_object
        for filter in reversed(self._filters):
            current_data = filter.get_filtered_object(current_data)
        return current_data
