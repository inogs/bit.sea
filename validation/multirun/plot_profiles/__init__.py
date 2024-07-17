from collections import namedtuple
from os import path
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np


from commons.mask import Mask
from commons.submask import SubMask
from commons.Timelist import TimeList
from commons import timerequestors
from commons import season
from utilities.mpi_serial_interface import get_mpi_communicator
from .tools.read_config import Config, InvalidConfigFile, PlotConfig, \
    DataDirSource, DepthProfilesOptions, TimeSeriesOptions, OutputOptions
from .tools.depth_profile_algorithms import get_depth_profile_plot_grid, \
    DepthProfileAlgorithm

try:
    from math import gcd, lcm
except ImportError:
    from math import gcd

    def lcm(a, b):
        return abs(a * b) // gcd(a, b)


# A BasinSliceVolume describes a slices among vertical levels of a specific
# basin, i.e., for example, all the cells between 10m and 30m in the Adriatic
# Sea (accordingly to a specific meshmask). The `levels` field should be a
# tuple with two elements: the start and the end of the interval in meters.
BasinSlice = namedtuple('BasinSlice', ['basin_index', 'meshmask', 'levels'])


class PlotDrawer:
    def __init__(self, config: Config, variable, meshmask_objects=None,
                 average_volume_weights=None):
        self._config = config
        self._plots = tuple(config.plots)
        self._levels = tuple(config.time_series_options.levels)
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
            if self._variable not in plot.variables:
                continue
            self._loaded_data[plot.name] = plot.get_plot_data(self._variable)

    def _get_plot_mask(self, plot):
        real_path = path.realpath(plot.source.meshmask)
        if real_path in self._meshmask_objects:
            return self._meshmask_objects[real_path]
        return Mask(
            plot.source.meshmask,
            maskvarname=plot.source.mask_var_name,
            loadtmask=True
        )

    def _plot_time_series(self, axis_dict, basin_index: int, basin,
                          is_2d: bool):
        elements_in_legend = False

        levels = self._levels
        if is_2d:
            levels = (0,)

        for plot in self._plots:
            if not plot.draw_time_series:
                continue

            if self._variable not in plot.variables:
                continue

            plot_meshmask = self._get_plot_mask(plot)

            if plot.name in self._loaded_data:
                plot_data = self._loaded_data[plot.name]
            else:
                plot_data = plot.get_plot_data(self._variable)

            for i, level in enumerate(levels):
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
                        warn(
                            'Average weights for level slice {} have not been '
                            'precomputed'.format(requested_slice)
                        )
                        submask = SubMask(basin, maskobject=plot_meshmask)
                        average_weights = compute_slice_volume(
                            level_slice,
                            plot_meshmask,
                            submask
                        )

                    if np.sum(np.abs(average_weights)) < 1e-12:
                        continue

                    with plot_data:
                        y_data_domain = plot_data.get_values(
                            time_steps=slice(None),
                            basin=basin_index,
                            level_index=level_slice
                        )
                    plot_y_data = np.average(
                        y_data_domain,
                        axis=-1,
                        weights=average_weights
                    )
                else:
                    level_index = plot_meshmask.getDepthIndex(level)
                    with plot_data:
                        plot_y_data = plot_data.get_values(
                            time_steps=slice(None),
                            basin=basin_index,
                            level_index=level_index
                        )
                current_axis = axis_dict['L{}'.format(i)]
                plt.sca(current_axis)

                # x and y labels
                if self._variable in self._config.variable_labels:
                    var_label = self._config.variable_labels[self._variable]
                    use_latex = False
                    if var_label.startswith('LaTeX:'):
                        use_latex = True
                        var_label = var_label[len('LaTeX:'):]
                    plt.ylabel(var_label, usetex=use_latex)

                # If we are drawing the last plot, add also the x_label
                if i == len(levels) - 1:
                    if self._config.time_series_options.x_label is not None:
                        plt.xlabel(self._config.time_series_options.x_label)

                # If there are no data, we skip this plot
                if np.all(np.isnan(plot_y_data)):
                    continue

                plot_kwargs = {}
                for field in ('alpha', 'color', 'zorder'):
                    if getattr(plot, field) is not None:
                        plot_kwargs[field] = getattr(plot, field)
                if plot.alpha_for_time_series is not None:
                    plot_kwargs['alpha'] = plot.alpha_for_time_series

                if plot.legend is not None:
                    plot_kwargs['label'] = plot.legend
                    elements_in_legend = True

                current_axis.plot(
                    plot_x_data,
                    plot_y_data,
                    **plot_kwargs
                )

        # Now we add the legends to each plot
        show_legend_flag = self._config.time_series_options.show_legend
        if show_legend_flag != "no" and elements_in_legend:
            only_one_legend = show_legend_flag in ("top", "bottom")

            axis_iterable = [
                axis_dict['L{}'.format(i)] for i in range(len(levels))
            ]

            if show_legend_flag == "bottom":
                axis_iterable = reversed(axis_iterable)

            for current_axis in axis_iterable:
                plt.sca(current_axis)
                if len(current_axis.lines) > 0:
                    plt.legend()
                    if only_one_legend:
                        break

        # Add the title to each plot; we need to do this here so the limits
        # of the plots are already defined
        for i, level in enumerate(levels):
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

            if self._config.time_series_options.show_depth:
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

        # Fix the x_ticks of the last axis
        last_l_axis = axis_dict['L{}'.format(len(levels) - 1)]
        if self._config.time_series_options.x_ticks_rotation is not None:
            for label in last_l_axis.get_xticklabels(which='major'):
                label.set_ha('right')
                label.set_rotation(
                    self._config.time_series_options.x_ticks_rotation
                )

    def _plot_depth_profile(self, axis_dict, basin_index: int):
        current_axis = axis_dict['P_0_0']

        elements_in_legend = False

        ytick_labels = self._config.depth_profiles_options.depth_ticks
        min_y = self._config.depth_profiles_options.min_depth
        max_y = self._config.depth_profiles_options.max_depth
        current_axis.set_yticks(ytick_labels)
        current_axis.set_ylim([min_y, max_y])
        current_axis.grid()
        current_axis.invert_yaxis()

        for plot in self._plots:
            if not plot.draw_depth_profile:
                continue
            if self._variable not in plot.variables:
                continue

            plot_meshmask = self._get_plot_mask(plot)

            if plot.name in self._loaded_data:
                plot_data = self._loaded_data[plot.name]
            else:
                plot_data = plot.get_plot_data(self._variable)

            with plot_data:
                plot_x_data = np.ma.mean(
                    plot_data.get_values(
                        time_steps=slice(None),
                        basin=basin_index,
                        level_index=slice(None)
                    ),
                    axis=0
                )
            plot_y_data = plot_meshmask.zlevels

            plot_kwargs = {}
            for field in ('alpha', 'color', 'zorder'):
                if getattr(plot, field) is not None:
                    plot_kwargs[field] = getattr(plot, field)

            if plot.legend is not None:
                plot_kwargs['label'] = plot.legend
                elements_in_legend = True

            current_axis.plot(
                plot_x_data,
                plot_y_data,
                **plot_kwargs
            )

        show_legend_flag = self._config.depth_profiles_options.show_legend
        if show_legend_flag != "no" and elements_in_legend:
            current_axis.legend()

        # Fix the x_ticks
        if self._config.depth_profiles_options.x_ticks_rotation is not None:
            for label in current_axis.get_xticklabels(which='major'):
                label.set_ha('right')
                label.set_rotation(
                    self._config.depth_profiles_options.x_ticks_rotation
                )

    def _plot_seasonal_depth_profile(self, axis_dict, basin_index: int):
        season_obj = season.season()
        elements_in_legend = False
        # x and y labels
        var_label = None
        use_latex = False
        if self._variable in self._config.variable_labels:
            var_label = self._config.variable_labels[self._variable]
            if var_label.startswith('LaTeX:'):
                use_latex = True
                var_label = var_label[len('LaTeX:'):]

        # Here we save which axes are on the last row
        bottom_axes = []

        for season_ind, season_str in enumerate(season_obj.SEASON_LIST_NAME):
            season_req = timerequestors.Clim_season(season_ind, season_obj)

            d_mode = self._config.depth_profiles_options.mode.config['mode']
            if d_mode == 'square':
                pi = season_ind // 2
                pj = season_ind % 2
            elif d_mode == 'inline':
                pi = 0
                pj = season_ind
            else:
                raise ValueError(
                    'Invalid display mode: "{}"'.format(d_mode)
                )

            current_axis = axis_dict[f'P_{pi}_{pj}']

            if d_mode == 'inline' or d_mode == 'square' and pi == 1:
                bottom_axes.append(current_axis)

            ytick_labels = self._config.depth_profiles_options.depth_ticks
            min_y = self._config.depth_profiles_options.min_depth
            max_y = self._config.depth_profiles_options.max_depth
            current_axis.set_ylim([min_y, max_y])
            current_axis.set_yticks(ytick_labels)
            current_axis.grid()
            current_axis.invert_yaxis()
            current_axis.set_title(season_str)

            if d_mode == 'square':
                if pi == 0:
                    current_axis.xaxis.set_tick_params(
                        which='both',
                        labelbottom=False,
                        labeltop=False
                    )
                if pi == 1:
                    if var_label is not None:
                        current_axis.set_xlabel(var_label, usetex=use_latex)

            if d_mode == 'inline':
                if pj > 0:
                    current_axis.set_yticklabels([])
                if var_label is not None:
                    current_axis.set_xlabel(var_label, usetex=use_latex)

            for plot in self._plots:

                if not plot.draw_depth_profile:
                    continue
                if self._variable not in plot.variables:
                    continue

                plot_meshmask = self._get_plot_mask(plot)

                if plot.name in self._loaded_data:
                    plot_data = self._loaded_data[plot.name]
                else:
                    plot_data = plot.get_plot_data(self._variable)

                with plot_data:

                    time_list = TimeList(plot_data.get_time_steps())

                    s_ind, s_weights = time_list.select(season_req)

                    plot_x_all_data = plot_data.get_values(
                        time_steps=slice(None), basin=basin_index,
                        level_index=slice(None))

                    plot_x_season_data = np.sum(
                        plot_x_all_data[s_ind, :] * s_weights.reshape((-1, 1)),
                        axis=0
                    ) / np.sum(s_weights)

                plot_y_data = plot_meshmask.zlevels

                plot_kwargs = {}
                for field in ('alpha', 'color', 'zorder'):
                    if getattr(plot, field) is not None:
                        plot_kwargs[field] = getattr(plot, field)

                if plot.legend is not None:
                    plot_kwargs['label'] = plot.legend
                    elements_in_legend = True

                current_axis.plot(
                    plot_x_season_data,
                    plot_y_data,
                    **plot_kwargs
                )

            show_legend_flag = self._config.depth_profiles_options.show_legend
            if show_legend_flag != "no" and elements_in_legend:
                current_axis.legend()

        # Fix the x_ticks
        if self._config.depth_profiles_options.x_ticks_rotation is not None:
            for axis in bottom_axes:
                for label in axis.get_xticklabels(which='major'):
                    label.set_ha('right')
                    label.set_rotation(
                        self._config.depth_profiles_options.x_ticks_rotation
                    )

    def plot(self, basin_index, basin, **fig_kw):
        if self.is_empty():
            raise ValueError('Required a plot from an empty PlotDrawer')

        draw_time_series = self._draw_time_series
        draw_depth_profile = self._draw_depth_profile
        levels = self._levels

        # We check if this is a 2d or a 3d plot
        is_2d = False
        if len(self._loaded_data) != 0:
            first_plot = list(self._loaded_data.values())[0]
            is_2d = first_plot.is_2d()
        else:
            for plot in self._plots:
                if self._variable not in plot.variables:
                    continue
                is_2d = plot.get_plot_data(self._variable).is_2d()
                break

        if is_2d:
            draw_time_series = True
            draw_depth_profile = False
            levels = (0, )

        # Building the grid of the plots
        # In this case, the plot grid is made only from the depth profile plots
        if not draw_time_series:
            plot_grid_rows, plot_grid_columns = get_depth_profile_plot_grid(
                self._config.depth_profiles_options.mode
            )
            if plot_grid_rows == 1 and plot_grid_columns == 1:
                fig, axis = plt.subplot(**fig_kw)
                axis_dict = {'P_0_0': axis}
            else:
                plot_structure = []
                for row in range(plot_grid_rows):
                    plot_structure.append(
                        ['P_{}_{}'.format(row, col) for col in
                         range(plot_grid_columns)]
                    )
                fig, axis_dict = plt.subplot_mosaic(
                    plot_structure,
                    **fig_kw
                )
        # Now the more general case
        else:
            if draw_depth_profile:
                dp_grid_rows, dp_grid_columns = get_depth_profile_plot_grid(
                    self._config.depth_profiles_options.mode
                )
                ts_ratio, dp_ratio = self._config.output_options.fig_ratio
            else:
                dp_grid_rows, dp_grid_columns = 1, 1
                ts_ratio, dp_ratio = 1, 1
            grid_rows = lcm(len(levels), dp_grid_rows)
            column_unit = dp_grid_columns // gcd(dp_ratio, dp_grid_columns)
            ts_cols_per_plot = ts_ratio * column_unit
            dp_cols_per_plot = dp_ratio * column_unit // dp_grid_columns

            plot_structure = []
            for row in range(grid_rows):
                current_level = row // (grid_rows // len(levels))
                current_dp_row = row // (grid_rows // dp_grid_rows)
                current_row = ['L{}'.format(current_level)] * ts_cols_per_plot
                if draw_depth_profile:
                    current_row.extend(
                        ['P_{}_{}'.format(current_dp_row, c)
                         for c in range(dp_grid_columns)
                         for _ in range(dp_cols_per_plot)]
                    )
                plot_structure.append(current_row)

            fig, axis_dict = plt.subplot_mosaic(
                plot_structure,
                **fig_kw
            )

            # Share the x-axis between each L plot and the last one
            last_l_axis = axis_dict['L{}'.format(len(levels) - 1)]
            for i in range(len(levels[:-1])):
                current_axis = axis_dict['L{}'.format(i)]
                current_axis.sharex(last_l_axis)
                current_axis.xaxis.set_tick_params(
                    which='both',
                    labelbottom=False,
                    labeltop=False
                )

        # If we have more than one depth profiles, we share all their axes
        if draw_depth_profile:
            dp_grid_rows, dp_grid_columns = get_depth_profile_plot_grid(
                self._config.depth_profiles_options.mode
            )
            if dp_grid_rows != 1 or dp_grid_columns != 1:
                first_axis = axis_dict['P_0_0']
                for pi in range(dp_grid_rows):
                    for pj in range(dp_grid_columns):
                        if pi == 0 and pj == 0:
                            continue
                        current_axis = axis_dict['P_{}_{}'.format(pi, pj)]
                        current_axis.sharex(first_axis)
                        current_axis.sharey(first_axis)

        if draw_time_series:
            self._plot_time_series(axis_dict, basin_index, basin, is_2d)
        if draw_depth_profile:
            if self._config.depth_profiles_options.mode.algorithm is \
                    DepthProfileAlgorithm.SEASONAL:
                self._plot_seasonal_depth_profile(axis_dict, basin_index)
            else:
                self._plot_depth_profile(axis_dict, basin_index)
        return fig


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

    n_cells = int(np.sum(basin_submask.mask[level_slice, :]))
    # If there are no cells inside this basin for the current slice, then there
    # is no need to compute its size
    if n_cells == 0:
        return np.zeros((e3t.shape[0],), dtype=np.float32)

    slice_volume = np.sum(
        e1t * e2t * e3t,
        axis=(1, 2),
        where=basin_submask.mask[level_slice, :]
    )

    return slice_volume


def draw_profile_plots(config: Config, basins, output_dir=None):
    # If no output_dir has been received, we use the output_dir that is saved
    # inside the config object
    if output_dir is None:
        output_dir = config.output_options.output_dir

    if output_dir is None:
        raise InvalidConfigFile(
            'output_dir has not been specified in the output section of the '
            'config file. Specify an output directory in the config file or '
            'submit the output_dir argument to this function.'
        )

    # Check if at least one of the specified levels requires to compute an
    # average
    averages_in_levels = False
    for level in config.time_series_options.levels:
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
                maskvarname=plot.source.mask_var_name,
                loadtmask=averages_in_levels
            )

    # And now we prepare the weights for the levels averages
    average_volume_weights = {}
    if averages_in_levels:
        for meshmask_path, meshmask_object in meshmask_objects.items():
            for basin_index, basin in enumerate(basins):
                submask = SubMask(basin, maskobject=meshmask_object)
                for level in config.time_series_options.levels:
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
    communicator = get_mpi_communicator()
    rank = communicator.Get_rank()
    comm_size = communicator.size
    for var_name in variables[rank::comm_size]:
        plot_drawer = PlotDrawer(
            config,
            var_name,
            meshmask_objects,
            average_volume_weights
        )
        if plot_drawer.is_empty():
            continue
        plot_drawer.load_data()

        for basin_index, basin in enumerate(basins):
            outfile_name = config.output_options.output_name.replace(
                '${VAR}',
                var_name
            )
            outfile_name = outfile_name.replace(
                '${BASIN}',
                basin.name
            )
            outfile_path = output_dir / outfile_name

            fig = plot_drawer.plot(
                basin_index=basin_index,
                basin=basin,
                dpi=config.output_options.dpi,
                figsize=config.output_options.fig_size
            )
            fig.suptitle('{} {}'.format(var_name, basin.name))

            fig.savefig(outfile_path)
            plt.close(fig)


__all__ = [
    PlotConfig, DataDirSource, DepthProfilesOptions, TimeSeriesOptions,
    OutputOptions, Config, draw_profile_plots
]
