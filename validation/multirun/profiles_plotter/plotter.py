from argparse import ArgumentParser
from os import path
from collections import namedtuple
from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np

from basins import V2
from commons.mask import Mask
from commons.submask import SubMask
from tools.read_config import read_config_from_file

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.size
    isParallel = True
except ModuleNotFoundError:
    rank = 0
    nranks = 1
    isParallel = False


MAIN_DIR = path.dirname(path.realpath(__file__))
CONFIG_FILE = path.join(MAIN_DIR, 'config.yaml')


# A BasinSliceVolume describes a slices among vertical levels of a specific
# basin, i.e., for example, all the cells between 10m and 30m in the Adriatic
# Sea (accordingly to a specific meshmask). The levels field should be a tuple
# with two elements: the start and the end of the interval in meters.
BasinSlice = namedtuple('BasinSlice', ['basin_index', 'meshmask', 'levels'])


class PlotDrawer:
    def __init__(self, plots: Iterable, variable, levels: Iterable,
                 meshmask_objects=None, average_volume_weights=None):
        self._plots = tuple(plots)
        self._levels = tuple(levels)
        self._variable = variable

        if meshmask_objects is None:
            self._meshmask_objects = {}
        else:
            self._meshmask_objects = meshmask_objects.copy()
        if average_volume_weights is None:
            self._average_volume_weights = {}
        else:
            self._average_volume_weights = average_volume_weights

        draw_depth_profile = False
        draw_time_series = False
        for plot in self._plots:
            if variable not in plot.variables:
                continue
            if plot.draw_depth_profile:
                draw_depth_profile = True
            if plot.draw_time_series:
                draw_time_series = True

        self._draw_depth_profile = draw_depth_profile
        self._draw_time_series = draw_time_series

        if len(self._levels) == 0:
            self._draw_time_series = False

        self._empty = not (self._draw_depth_profile or self._draw_time_series)
        self._loaded_data = {}

    def is_empty(self):
        return self._empty

    @property
    def levels(self):
        return self._levels

    def load_data(self):
        self._loaded_data = {}
        for plot in self._plots:
            self._loaded_data[plot.name] = plot.get_plot_data(self._variable)

    def _get_plot_mask(self, plot):
        real_path = path.realpath(plot.source.meshmask)
        if real_path in self._meshmask_objects:
            return self._meshmask_objects[real_path]
        return Mask(plot.source.meshmask, loadtmask=True)

    def _plot_time_series(self, axis_dict, basin_index: int, basin,
                          indicator=0):
        for plot in self._plots:
            if not plot.draw_time_series:
                continue

            plot_meshmask = self._get_plot_mask(plot)

            if plot.name in self._loaded_data:
                plot_data = self._loaded_data[plot.name]
            else:
                plot_data = plot.get_plot_data(self._variable)

            for i, level in enumerate(self._levels):
                plot_x_data = plot_data.get_time_steps()

                # If level is a tuple, we have to compute the average between
                # its first number and the second
                if isinstance(level, tuple):
                    # Let us check if we already have the weights for this
                    # average
                    requested_slice = BasinSlice(
                        basin_index,
                        path.realpath(plot.source.meshmask),
                        level
                    )

                    level_slice = from_interval_to_slice(
                        level,
                        plot_meshmask
                    )

                    if requested_slice in self._average_volume_weights:
                        average_weights = \
                            self._average_volume_weights[requested_slice]
                    else:
                        submask = SubMask(basin, maskobject=plot_meshmask)
                        average_weights = compute_slice_volume(
                            level_slice,
                            plot_meshmask,
                            submask
                        )

                    with plot_data:
                        y_data_domain = plot_data.get_values(
                            time_steps=slice(None),
                            basin=basin_index,
                            level_index=level_slice,
                            indicator=indicator  # Zero is the mean
                        )
                    plot_y_data = np.average(
                        y_data_domain,
                        axis=plot_data.get_axis(
                            'level',
                            while_fixing={'basin', 'coasts'}
                        ),
                        weights=average_weights
                    )
                else:
                    level_index = plot_meshmask.getDepthIndex(level)
                    with plot_data:
                        plot_y_data = plot_data.get_values(
                            time_steps=slice(None),
                            basin=basin_index,
                            level_index=level_index,
                            indicator=indicator  # Zero is the mean
                        )

                if np.all(np.isnan(plot_y_data)):
                    continue

                plot_kwargs = {}
                for field in ('alpha', 'color', 'zorder'):
                    if getattr(plot, field) is not None:
                        plot_kwargs[field] = getattr(plot, field)
                if plot.alpha_for_time_series is not None:
                    plot_kwargs['alpha'] = plot.alpha_for_time_series

                current_axis = axis_dict['L{}'.format(i)]
                current_axis.plot(
                    plot_x_data,
                    plot_y_data,
                    label=plot.legend,
                    **plot_kwargs
                )

        # Add the title to each plot; we need to do this here so the limits
        # of the plots are already defined
        for i, level in enumerate(self._levels):
            current_axis = axis_dict['L{}'.format(i)]
            plot_left, plot_right = current_axis.get_xlim()
            plot_bottom, plot_top = current_axis.get_ylim()

            axis_title = 'Depth: '
            if isinstance(level, tuple):
                if level[0] is None:
                    start_average_str = '0m'
                else:
                    start_average_str = '{}m'.format(level[0])
                if level[-1] is None:
                    end_average_str = 'bottom'
                else:
                    end_average_str = '{}m'.format(level[-1])
                depth_str = 'from {} to {}'.format(
                    start_average_str,
                    end_average_str
                )
            else:
                depth_str = '{}m'.format(level)
            axis_title += depth_str

            title_x_pos = plot_left + (plot_right - plot_left) * 0.03
            title_y_pos = plot_top - (plot_top - plot_bottom) * 0.05
            current_axis.text(
                title_x_pos,
                title_y_pos,
                axis_title,
                fontsize=10,
                ha='left',
                va='top'
            )

    def _plot_depth_profile(self, axis_dict, basin_index: int):
        current_axis = axis_dict['P']

        ytick_labels = (0, 200, 400, 600, 800, 1000)
        current_axis.set_ylim([0, 1000])
        current_axis.set_yticks(ytick_labels)
        current_axis.grid()
        current_axis.invert_yaxis()

        for plot in self._plots:
            if not plot.draw_depth_profile:
                continue

            plot_meshmask = self._get_plot_mask(plot)

            if plot.name in self._loaded_data:
                plot_data = self._loaded_data[plot.name]
            else:
                plot_data = plot.get_plot_data(self._variable)

            with plot_data:
                plot_x_data = np.mean(
                    plot_data.get_values(
                        time_steps=slice(None),
                        basin=basin_index,
                        level_index=slice(None),
                        indicator=0
                    ),
                    axis=0
                )
            plot_y_data = plot_meshmask.zlevels

            plot_kwargs = {}
            for field in ('alpha', 'color', 'zorder'):
                if getattr(plot, field) is not None:
                    plot_kwargs[field] = getattr(plot, field)

            current_axis.plot(
                plot_x_data,
                plot_y_data,
                label=plot.legend,
                **plot_kwargs
            )
        current_axis.legend()

    def plot(self, basin_index, basin, **fig_kw):
        if self.is_empty():
            raise ValueError('Required a plot from an empty PlotDrawer')

        if not self._draw_time_series:
            fig, axis = plt.subplot(**fig_kw)
            axis_dict = {'P': axis}
        else:
            plot_structure = [
                ['L{}'.format(i)] for i in range(len(self._levels))
            ]
            if self._draw_depth_profile:
                for plot_row in plot_structure:
                    plot_row.append('P')
            fig, axis_dict = plt.subplot_mosaic(
                plot_structure,
                **fig_kw
            )

            # Share the x-axis between each L plot and the last one
            last_l_axis = axis_dict['L{}'.format(len(self._levels) - 1)]
            for i in range(len(self._levels[:-1])):
                current_axis = axis_dict['L{}'.format(i)]
                current_axis.sharex(last_l_axis)
                current_axis.xaxis.set_tick_params(
                    which='both',
                    labelbottom=False,
                    labeltop=False
                )
                # If there are more than 3 time series, hide all the x-ticks
                # on each plot beside the bottom one
                if len(self._levels[:-1]) > 3:
                    current_axis.xaxis.offsetText.set_visible(False)

        if self._draw_time_series:
            self._plot_time_series(axis_dict, basin_index, basin)
        if self._draw_depth_profile:
            self._plot_depth_profile(axis_dict, basin_index)

        return fig


def configure_argparse():
    parser = ArgumentParser()
    parser.description = (
        "The MultirunProfilePlotter is a script to plot temporal series and "
        "depth profiles of several different runs of one or more models. "
        "You may configure the behaviour of the script by setting the values "
        "of the config.yaml file in the main directory of this script or "
        "writing another configuration file yourself."
    )
    parser.add_argument(
        'config_file',
        type=str,
        nargs='?',
        default=CONFIG_FILE,
        help='The path of the config file used by this script. By default, it '
             'uses {}'.format(CONFIG_FILE)
    )
    return parser.parse_args()


def from_interval_to_slice(level_interval: tuple, meshmask: Mask):
    """
    Given an interval in meters, return the slice that identifies that interval
    as z-levels of the model.

    :param level_interval: a tuple with two numbers (in meters)
    :param meshmask: the meshmask of the model
    :return: a slice that, when applied to the levels of the model, returns the
    cells whose depth is among the values specified in level_interval
    """
    if len(level_interval) != 2:
        raise ValueError(
            'Levels must be a number of a tuple with two '
            'elements (the start and the end of the interval '
            'where an average will be computed). Submitted a '
            'tuple with {} elements: {}'.format(
                len(level_interval),
                level_interval
            )
        )
    if level_interval[0] is None:
        start_level_index = None
    else:
        start_level_index = meshmask.getDepthIndex(
            level_interval[0]
        )
    if level_interval[-1] is None:
        end_level_index = None
    else:
        end_level_index = meshmask.getDepthIndex(level_interval[-1])
        # This is because we want to include the last level
        # inside the average
        end_level_index += 1

    level_slice = slice(start_level_index, end_level_index)
    return level_slice


def compute_slice_volume(level_slice: slice, meshmask: Mask,
                         basin_submask: SubMask) -> np.ndarray:
    """
    Compute the volume of each level for a specific mask and basin

    :param level_slice: Level slice is a slice that selects the levels where we
    want to compute the volumes. This is a slice for the levels, it is not in
    meters
    :param meshmask: The meshmask of the model
    :param basin_submask: The submask of the basin
    :return:
    A 1D array `v` such that `v[i]` is the volume of the i-th level of the
    slice
    """
    e1t = meshmask.e1t
    e2t = meshmask.e2t
    e3t = meshmask.e3t[level_slice, :]

    slice_volume = np.sum(
        e1t * e2t * e3t,
        axis=(1, 2),
        where=basin_submask.mask[level_slice, :]
    )
    return slice_volume


def main():
    args = configure_argparse()

    config = read_config_from_file(args.config_file)

    # Check if at least one of the specified levels requires to compute an
    # average
    averages_in_levels = False
    for level in config.levels:
        if isinstance(level, tuple):
            averages_in_levels = True
            break

    # Create a list of all the variables (for all the plots) and of all the
    # meshmasks
    variables = set()
    meshmask_objects = {}
    for plot in config.plots:
        variables.update(plot.variables)
        meshmask_path = path.realpath(plot.source.meshmask)
        # Here we read in advance the meshmask objects, so that we can share
        # the information among different plots (if they use the same
        # meshmask). If we need to compute the volumes of the levels
        # (because there are some averages that must be computed), we also read
        # the t-mask
        if meshmask_path not in meshmask_objects:
            meshmask_objects[meshmask_path] = Mask(
                meshmask_path,
                loadtmask=averages_in_levels
            )

    # And now we prepare the weights for the levels averages
    average_volume_weights = {}
    if averages_in_levels:
        for meshmask_path, meshmask_object in meshmask_objects.items():
            for basin_index, basin in enumerate(V2.P):
                submask = SubMask(basin, maskobject=meshmask_object)
                for level in config.levels:
                    if not isinstance(level, tuple):
                        continue
                    level_slice = from_interval_to_slice(
                        level,
                        meshmask_object
                    )
                    selection_weights = compute_slice_volume(
                        level_slice,
                        meshmask_object,
                        submask
                    )
                    basin_slice = BasinSlice(basin_index, meshmask_path, level)
                    average_volume_weights[basin_slice] = selection_weights

    # Now we sort the variables, so we ensure that the list is in the same
    # order on all the processes
    variables = tuple(sorted(list(variables)))

    # Here we produce the plots
    for var_name in variables[rank::nranks]:
        plot_drawer = PlotDrawer(
            config.plots,
            var_name,
            config.levels,
            meshmask_objects,
            average_volume_weights
        )
        if plot_drawer.is_empty():
            continue
        plot_drawer.load_data()

        for basin_index, basin in enumerate(V2.P):
            outfile_name = 'Multirun_Profiles.{}.{}.png'.format(
                var_name,
                basin.name
            )
            outfile_path = config.output_dir / outfile_name

            fig = plot_drawer.plot(
                basin_index=basin_index,
                basin=basin,
                figsize=(10, 10)
            )
            fig.suptitle('{} {}'.format(var_name, basin.name))

            fig.savefig(outfile_path)
            plt.close(fig)


if __name__ == '__main__':
    main()
