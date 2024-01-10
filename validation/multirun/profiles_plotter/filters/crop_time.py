from datetime import datetime
import re
from typing import Union

import numpy as np

from filters.filter import SimpleFilter, FilteredObject
from tools.data_object import DataObject


class AllPointsCropped(ValueError):
    pass


class TimeCroppedObject(FilteredObject):
    def __init__(self, data_object: DataObject,
                 start_time: Union[datetime, None] = None,
                 end_time: Union[datetime, None] = None):
        super().__init__(data_object)

        old_time_steps = data_object.get_time_steps()

        start_index = 0
        if start_time is not None:
            while start_index < len(old_time_steps) and \
                    old_time_steps[start_index] < start_time:
                start_index += 1
            if old_time_steps[start_index] < start_time:
                raise AllPointsCropped(
                    'No points survived after the cropping: they are all '
                    'before {}'.format(start_time)
                )
        start_time = old_time_steps[start_index]

        end_index = len(old_time_steps)
        if end_time is not None:
            while end_index > 1 and old_time_steps[end_index - 1] > end_time:
                end_index -= 1
            if old_time_steps[end_index - 1] > end_time:
                raise AllPointsCropped(
                    'No points survived after the cropping: they are all '
                    'after {}'.format(end_time)
                )
        end_time = old_time_steps[end_index - 1]

        self._start_time, self._end_time = start_time, end_time
        self._start_index, self._end_index = start_index, end_index

        self._time_steps = old_time_steps[start_index:end_index]

    def get_time_steps(self):
        return self._time_steps

    @property
    def start_time(self):
        return self._start_time

    @property
    def end_time(self):
        return self._end_time

    def get_values(self, time_steps, basin, level_index, indicator, coasts=0):
        if isinstance(time_steps, slice):
            if time_steps.start is None:
                original_slice_start = self._start_index
            else:
                time_steps_start = time_steps.start
                if time_steps_start < 0:
                    time_steps_start += len(self._time_steps)
                if time_steps_start < 0:
                    original_slice_start = self._start_index
                else:
                    original_slice_start = self._start_index + time_steps_start

            if time_steps.stop is None:
                original_slice_stop = self._end_index
            else:
                time_steps_stop = time_steps.stop
                if time_steps_stop < 0:
                    time_steps_stop += len(self._time_steps)
                if time_steps_stop < 0:
                    original_slice_stop = 0
                else:
                    original_slice_stop = self._start_index + time_steps_stop

            original_time_slice = slice(
                original_slice_start,
                original_slice_stop,
                time_steps.step
            )
        else:
            # Here we can suppose that time_steps is an integer
            if time_steps >= len(self._time_steps) or \
                    time_steps < - len(self._time_steps):
                raise IndexError(
                    'Index {} is out of bounds for time axis with size '
                    '{}'.format(time_steps, len(self._time_steps))
                )
            if time_steps < 0:
                time_steps += len(self._time_steps)
            original_time_slice = time_steps + self._start_index

        with self._original_data:
            original_data = self._original_data.get_values(
                time_steps=original_time_slice,
                basin=basin,
                level_index=level_index,
                indicator=indicator,
                coasts=coasts
            )

        return original_data


class TimeCroppingFilter(SimpleFilter):
    TIME_FORMAT = '%Y%m%d-%H:%M:%S'

    def __init__(self, start_time: Union[datetime, None] = None,
                 end_time: Union[datetime, None] = None):
        self._start_time = start_time
        self._end_time = end_time

    def get_filtered_object(self, data_object) -> TimeCroppedObject:
        return TimeCroppedObject(
            data_object=data_object,
            start_time=self._start_time,
            end_time=self._end_time
        )

    @staticmethod
    def initialize_from_string(args_str):
        error_message = \
            'A {} requires two arguments: a starting date and an ending ' \
            'date, in the following format: {}.'.format(
                TimeCroppingFilter.__name__,
                TimeCroppingFilter.TIME_FORMAT
            )
        if args_str is None:
            raise ValueError(error_message + ' None received.')

        args_split = args_str.split(',')
        if len(args_split) < 2:
            raise ValueError(error_message + 'Received only one argument.')
        if len(args_split) > 2:
            raise ValueError(error_message + 'Received {} arguments.'.format(
                len(args_split)
            ))
        args_split = tuple(p.strip() for p in args_split)

        start_time = datetime.strptime(
            args_split[0],
            TimeCroppingFilter.TIME_FORMAT
        )
        end_time = datetime.strptime(
            args_split[1],
            TimeCroppingFilter.TIME_FORMAT
        )
        return TimeCroppingFilter(start_time, end_time)
