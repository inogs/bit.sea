from datetime import timedelta
import re

import numpy as np

from filters.filter import SimpleFilter, FilteredObject
from plot_inputs import PlotInputData


class MovingAverageFilteredObject(FilteredObject):
    def __init__(self, original_input: PlotInputData, window_size: int):
        super().__init__(original_input)

        self._window_size = int(window_size)

        # How many points do we loose because of the average (the ones on the
        # boundary of the time series)
        lost_points = self._window_size - 1

        original_times = self._input_data.get_time_steps()

        if self._window_size % 2 == 1:
            lost_per_side = lost_points // 2
            self._time_steps = tuple(
                original_times[lost_per_side: - lost_per_side]
            )
        else:
            lost_per_side = (lost_points - 1) // 2
            new_time_steps = []
            averaged_times = original_times[lost_per_side:- lost_per_side]
            for i, t1, in enumerate(averaged_times[:-1]):
                t2 = averaged_times[i + 1]
                delta_time = t2 - t1
                new_time_steps.append(t1 + delta_time // 2)
            self._time_steps = tuple(new_time_steps)
            assert len(self._time_steps) == len(original_times) - lost_points

    def get_time_steps(self):
        return self._time_steps

    def get_values(self, time_steps, basin, level_index):
        fixed_axis = {'coasts', 'indicator'}
        if not isinstance(basin, slice):
            fixed_axis.add('basin')
        if not isinstance(level_index, slice):
            fixed_axis.add('level')

        time_axis = self._data_object.get_axis('time', fixed_axis)

        if isinstance(time_steps, slice):
            if time_steps.stop is None:
                new_slice_stop = None
            else:
                if time_steps.stop < 0:
                    new_slice_stop = time_steps.stop
                else:
                    new_slice_stop = time_steps.stop + self._window_size - 1
            time_slice = slice(time_steps.start, new_slice_stop)
        else:
            # Here we can suppose that time_steps is an integer
            if time_steps >= 0:
                time_slice_start = time_steps
                time_slice_stop = time_steps + self._window_size
            else:
                time_slice_start = time_steps - self._window_size
                time_slice_stop = time_steps
            time_slice = slice(time_slice_start, time_slice_stop)

        original_data = self._input_data.get_values(
            time_steps=time_slice,
            basin=basin,
            level_index=level_index
        )

        if isinstance(time_steps, int):
            if original_data.shape[time_axis] != self._window_size:
                raise IndexError(
                    'Index {} outside valid range'.format(time_steps)
                )

            return np.mean(original_data, axis=time_axis)

        all_slice_none = [
            slice(None) for _ in range(len(original_data.shape))
        ]

        time_index_geq_0 = all_slice_none.copy()
        time_index_geq_0[time_axis] = slice(1, None)
        time_index_geq_0 = tuple(time_index_geq_0)

        cut_first_time_indices = all_slice_none.copy()
        cut_first_time_indices[time_axis] = slice(self._window_size - 1, None)
        cut_first_time_indices = tuple(cut_first_time_indices)

        cut_last_time_indices = all_slice_none.copy()
        cut_last_time_indices[time_axis] = slice(None, - self._window_size)
        cut_last_time_indices = tuple(cut_last_time_indices)

        v1 = np.cumsum(
            original_data,
            axis=time_axis
        )[cut_first_time_indices]

        v2 = np.zeros_like(v1)
        v2[time_index_geq_0] = np.cumsum(
            original_data,
            axis=time_axis
        )[cut_last_time_indices]

        moving_average = (v1 - v2) / self._window_size

        return moving_average


class MovingAverageFilter(SimpleFilter):
    WINDOW_DESCRIPTION_MASK = re.compile(
        r'^(?P<size>\d+)\s{0,2}(?P<unit>([dDmMyY])?)\s*$'
    )

    def __init__(self, window_description: str):
        mask_match = self.WINDOW_DESCRIPTION_MASK.match(window_description)
        if mask_match is None:
            raise ValueError(
                'The string that describes the window size of the moving '
                'average filter must be an integer, optionally followed by '
                'one of the following letters: d, m or y'
            )
        self._window_size = int(mask_match.group('size'))
        if mask_match.group('unit') is not None:
            self._unit = mask_match.group('unit').lower()
        else:
            self._unit = None

    def get_filtered_object(self, input_data: PlotInputData):
        # If there are no units, we can simply return the filtered object as is
        if self._unit is None:
            return MovingAverageFilteredObject(
                input_data,
                window_size=self._window_size
            )

        # Otherwise, we need to estimate the new size of the window
        if self._unit == 'd':
            unit_time = timedelta(days=1)
        elif self._unit == 'm':
            unit_time = timedelta(days=30)
        else:
            unit_time = timedelta(days=365)

        time_steps = input_data.get_time_steps()
        time_index = 0
        for time_index, current_time in enumerate(time_steps):
            if current_time - time_steps[0] >= unit_time:
                break

        if time_index == 0:
            raise ValueError('No time steps found in the data_object')

        window_size = self._window_size * time_index
        return MovingAverageFilteredObject(input_data, window_size=window_size)

    @staticmethod
    def initialize_from_string(args_str):
        if args_str is None:
            raise ValueError(
                'A moving average filter requires a mandatory argument that '
                'must be the size of the average window'
            )

        args_split = args_str.split(',')
        if len(args_split) > 1:
            raise ValueError(
                'A moving average filter requires only one argument. '
                'Received: {}'.format(
                    ', '.join([k.strip() for k in args_split])
                )
            )
        window_desc = args_split[0].strip()
        return MovingAverageFilter(window_desc)
