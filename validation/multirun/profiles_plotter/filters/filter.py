from abc import ABC, abstractmethod
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

    @staticmethod
    @abstractmethod
    def initialize_from_string(args_str: Union[str, None]):
        raise NotImplementedError
