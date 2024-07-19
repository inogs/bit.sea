from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
import yaml
from math import gcd
from numbers import Real
from pathlib import Path
import re
from sys import version_info
from typing import Any, Dict, Literal, Mapping, Tuple, Union

from ..tools.data_object import DataObject, PickleDataObject
from ..plot_inputs import PlotInputData
from ..plot_inputs.single_line_plot import SingleLineInputData
from ..tools.depth_profile_algorithms import DEFAULT_DEPTH_PROFILE_MODE, \
    read_depth_profile_mode, DepthProfileMode
from ..filters.read_filter_description import read_filter_description

if version_info[1] < 9:
    from typing import Callable
else:
    from collections.abc import Callable


DEFAULT_OUTPUT_NAME = 'Multirun_Profiles.${VAR}.${BASIN}.png'
DEFAULT_DPI = 300
DEFAULT_FIG_SIZE = (10, 10)
DEFAULT_FIG_RATIO = (1, 1)

DEFAULT_COAST_INDEX = 1
DEFAULT_INDICATOR = 0

DEFAULT_MASK_VAR_NAME = 'tmask'

TICK_STR_MASK = re.compile(
    r'^\s*([Rr]ange|RANGE)\s*\((?P<arguments>.*)\)\s*$'
)


EXPECTED_FIELDS = {
    'sources',
    'variable_selections',
    'variable_labels',
    'plots',
    'levels',
    'time_series',
    'depth_profiles',
    'output'
}

OUTPUT_FIELDS = (
    'output_name',
    'output_dir',
    'dpi',
    'fig_size',
    'fig_ratio'
)


class InvalidConfigFile(Exception):
    pass


class InvalidPlotConfig(InvalidConfigFile):
    pass


@dataclass
class DataDirSource:
    path: Path
    meshmask: Path
    cost_index: Union[int, None] = None
    mask_var_name: str = DEFAULT_MASK_VAR_NAME


@dataclass
class OutputOptions:
    output_name: str = DEFAULT_OUTPUT_NAME
    output_dir: Union[Path, None] = None
    dpi: int = DEFAULT_DPI
    fig_size: Tuple[Real, Real] = DEFAULT_FIG_SIZE
    fig_ratio: Tuple[Real, Real] = DEFAULT_FIG_RATIO


@dataclass
class TimeSeriesOptions:
    levels: Tuple
    show_legend: Literal["no", "all", "bottom", "top"] = "no"
    show_depth: bool = True
    x_label: Union[str, None] = None
    x_ticks_rotation: Union[Real, None] = None
    show_grid: bool = False


@dataclass
class DepthProfilesOptions:
    mode: DepthProfileMode = DEFAULT_DEPTH_PROFILE_MODE
    show_legend: Literal["no", "all", "bottom", "top"] = "all"
    depth_ticks: Tuple[Real, ...] = tuple(range(0, 1001, 200))
    min_depth: Real = 0
    max_depth: Real = 1000
    x_ticks_rotation: Union[Real, None] = None
    show_y_ticks: Literal["all", "no", "left", "right"] = "all"
    y_ticks_position: Literal["left", "right"] = "left"


def read_number(number_str: str) -> Union[int, float]:
    try:
        return int(number_str)
    except ValueError:
        return float(number_str)


def read_output_dir(output_dir: str) -> Path:
    output_dir = Path(output_dir)

    # If output dir does not exist, create it
    if not output_dir.exists():
        if not output_dir.parent.exists():
            raise IOError('Directory {} does not exists'.format(output_dir))
        output_dir.mkdir(exist_ok=True)
    return output_dir


INTEGER_PAIR = re.compile(
    r'^\s*\(\s*(?P<first>\d+)\s*,\s*(?P<second>\d+)\s*\)\s*$'
)


def read_pair_of_integers(fig_size_str: str) -> Tuple[int, int]:
    str_match = INTEGER_PAIR.match(fig_size_str)
    if str_match is None:
        raise ValueError(
            'Invalid string for a pair of integers: {}'.format(fig_size_str)
        )
    first = int(str_match.group('first'))
    second = int(str_match.group('second'))

    return first, second


def read_time_series_options(option_dict: Dict[str, Any]) -> TimeSeriesOptions:
    show_legend_default = "no"
    show_depth_default = False
    x_label_default = None

    # Read levels, ensuring that they are a list of integers or floats
    if 'levels' not in option_dict:
        raise InvalidConfigFile(
            'No "levels" field found inside the "time_series" options'
        )

    levels_raw = option_dict['levels']
    levels = []
    if not isinstance(levels_raw, Iterable):
        raise InvalidConfigFile(
            'Field "levels" must be a list of levels (in meters)'
        )
    for raw_level in levels_raw:
        if isinstance(raw_level, str) and '-' in raw_level:
            level_split = raw_level.split('-')
            if len(level_split) > 2:
                raise InvalidConfigFile(
                    'Invalid level found: {}'.format(raw_level)
                )
            level_start_raw = level_split[0].strip()
            level_end_raw = level_split[1].strip()
            if level_start_raw == '':
                level_start = None
            else:
                level_start = read_number(level_start_raw)
            if level_end_raw == '':
                level_end = None
            else:
                level_end = read_number(level_end_raw)
            levels.append((level_start, level_end))
            continue
        levels.append(read_number(raw_level))
    levels = tuple(levels)

    for key in option_dict:
        if key not in ('levels', 'show_legend', 'show_depth', 'x_label'):
            raise ValueError(
                'Invalid key inside the time series option dictionary: '
                '{}'.format(key)
            )

    x_label = x_label_default
    if 'x_label' in option_dict:
        x_label = option_dict['x_label']

    show_legend = show_legend_default
    if 'show_legend' in option_dict:
        show_legend_raw = option_dict['show_legend']
        if show_legend_raw is None:
            show_legend = 'no'
        elif isinstance(show_legend_raw, bool):
            show_legend = 'all' if show_legend_raw else 'no'
        else:
            show_legend = str(show_legend_raw).strip().lower()
            if show_legend not in ('all', 'bottom', 'top', 'no'):
                raise ValueError(
                    'time_series:show_legend field must contain either "all" '
                    'or "bottom" or "top" or "no"; received {}'.format(
                        show_legend
                    )
                )

    show_depth = show_depth_default
    if 'show_depth' in option_dict:
        if not isinstance(option_dict['show_depth'], bool):
            if option_dict['show_depth'] is not None:
                raise ValueError(
                    'time_series:show_depth field must contain a boolean value'
                )
        if option_dict['show_depth'] is None:
            show_depth = False
        else:
            show_depth = bool(option_dict['show_depth'])

    return TimeSeriesOptions(
        levels=levels,
        show_legend=show_legend,
        show_depth=show_depth,
        x_label=x_label
    )


def read_ticks(ticks_raw: Union[None, str, list]) -> Tuple[Real, ...]:
    if ticks_raw is None:
        return ()
    if isinstance(ticks_raw, list):
        ticks_list = []
        for tick in ticks_raw:
            try:
                current_tick = read_number(tick)
            except ValueError:
                raise ValueError('Invalid tick value: {}'.format(tick))
            ticks_list.append(current_tick)
        return tuple(ticks_list)

    tick_str_match = TICK_STR_MASK.match(ticks_raw)
    if tick_str_match is None:
        raise ValueError('Unknown depth_ticks string: {}'.format(ticks_raw))

    ticks_args_raw = [
        t.strip() for t in tick_str_match.group('arguments').split(',')
    ]
    ticks_args = []
    for raw_tick in ticks_args_raw:
        try:
            ticks_args.append(read_number(raw_tick))
        except ValueError:
            raise ValueError(
                'Invalid arg inside range for depth_ticks: "{}"'.format(
                    raw_tick
                )
            )
    if len(ticks_args) == 1:
        start = 0
        end = ticks_args[0]
        step = 1
    elif len(ticks_args) == 2:
        start, end = ticks_args
        step = 1
    elif len(ticks_args) == 3:
        start, end, step = ticks_args
    else:
        raise ValueError(
            'A range can not accept {} arguments (inside depth_ticks '
            'field)'.format(len(ticks_args))
        )

    if step <= 0:
        raise ValueError('step must be a positive number')
    if end < start:
        raise ValueError(
            'The end value of the range must be bigger than the start'
        )

    ticks = []
    current_tick = start
    while current_tick < end:
        ticks.append(current_tick)
        current_tick += step
        if len(ticks) > 100:
            raise ValueError(
                'Value "{}" inside depth_ticks field generates too many'
                'ticks'.format(ticks_raw)
            )
    return tuple(ticks)


def read_depth_profiles_options(
            option_dict: Union[Dict[str, Any], None] = None
        ) -> DepthProfilesOptions:
    mode = DEFAULT_DEPTH_PROFILE_MODE
    show_legend_default = "all"
    depth_ticks_raw = 'range(0, 1001, 200)'

    if option_dict is None:
        depth_ticks = read_ticks(depth_ticks_raw)
        max_depth = depth_ticks[-1] if len(depth_ticks) > 0 else 1000
        min_depth = depth_ticks[0] if len(depth_ticks) > 0 else 0
        return DepthProfilesOptions(
            mode=mode,
            show_legend=show_legend_default,
            depth_ticks=read_ticks(depth_ticks_raw),
            min_depth=min_depth,
            max_depth=max_depth,
        )

    show_legend = show_legend_default
    if 'show_legend' in option_dict:
        show_legend_raw = option_dict['show_legend']
        if show_legend_raw is None:
            show_legend = 'no'
        elif isinstance(show_legend_raw, bool):
            show_legend = 'all' if show_legend_raw else 'no'
        else:
            show_legend = str(show_legend_raw).strip().lower()
            if show_legend not in ('all', 'yes', 'no'):
                raise ValueError(
                    'depth_profiles:show_legend field must contain either'
                    '"yes" or "no"; received {}'.format(
                        show_legend
                    )
                )
            if show_legend == 'yes':
                show_legend = "all"

    for key in option_dict:
        if key not in (
                'mode', 'show_legend', 'depth_ticks', 'min_depth', 'max_depth'
        ):
            raise ValueError(
                'Invalid key inside the depth profiles option dictionary: '
                '{}'.format(key)
            )
        if key == 'mode':
            mode = read_depth_profile_mode(option_dict[key])
        if key == 'depth_ticks':
            depth_ticks_raw = option_dict['depth_ticks']

    depth_ticks = read_ticks(depth_ticks_raw)
    max_depth = depth_ticks[-1] if len(depth_ticks) > 0 else 1000
    min_depth = depth_ticks[0] if len(depth_ticks) > 0 else 0

    if 'min_depth' in option_dict:
        try:
            min_depth = read_number(option_dict['min_depth'])
        except ValueError:
            raise ValueError(
                'Invalid value for "min_depth": {}'.format(
                    option_dict['min_depth']
                )
            )
    if 'max_depth' in option_dict:
        try:
            max_depth = read_number(option_dict['max_depth'])
        except ValueError:
            raise ValueError(
                'Invalid value for "max_depth": {}'.format(
                    option_dict['max_depth']
                )
            )

    return DepthProfilesOptions(
        mode=mode,
        show_legend=show_legend,
        depth_ticks=depth_ticks,
        min_depth=min_depth,
        max_depth=max_depth
    )


class PlotConfig:
    CONFIG_FIELDS = {
        'source',
        'variables',
        'coast_index',
        'indicator',
        'active',
        'color',
        'legend',
        'alpha',
        'alpha_for_time_series',
        'zorder',
        'filter',
        'draw_time_series',
        'draw_depth_profile'
    }

    PLOT_NAMES = set()

    def __init__(self, name: str, source: DataDirSource,
                 variables: Tuple[str, ...],
                 plot_builder: Callable[[DataObject], PlotInputData],
                 data_filter=None, active=True, color=None, alpha=None,
                 alpha_for_time_series=None, legend=None, zorder=None,
                 draw_time_series=True, draw_depth_profile=True):
        if name in self.PLOT_NAMES:
            raise ValueError(
                'A plot named "{}" has already been created'.format(name)
            )
        self.PLOT_NAMES.add(name)

        self.name: str = name
        self._source: DataDirSource = source
        self.variables: Tuple[str, ...] = variables
        self._plot_builder: Callable[[DataObject], PlotInputData] = \
            plot_builder

        self._active: bool = active
        self._filter = data_filter

        self._color = color
        self._alpha = alpha
        self._alpha_for_time_series = alpha_for_time_series
        self._legend: Union[str, None] = legend
        self._zorder: Union[int, None] = zorder

        self._draw_time_series: bool = bool(draw_time_series)
        self._draw_depth_profile: bool = bool(draw_depth_profile)

    def get_plot_data(self, variable):
        if variable not in self.variables:
            raise ValueError(
                'Plot "{}" does not contain variable "{}"'.format(
                    self.name,
                    variable
                )
            )

        data_object = PickleDataObject(self.source.path, variable)
        plot_input = self._plot_builder(data_object)
        if self._filter is None:
            return plot_input

        return self._filter.get_filtered_object(plot_input)

    @property
    def source(self) -> DataDirSource:
        return self._source

    def is_active(self) -> bool:
        return self._active

    @property
    def color(self):
        return self._color

    @property
    def alpha(self):
        return self._alpha

    @property
    def alpha_for_time_series(self):
        return self._alpha_for_time_series

    @property
    def legend(self) -> Union[str, None]:
        return self._legend

    @property
    def zorder(self) -> int:
        return self._zorder

    @property
    def draw_time_series(self) -> bool:
        return self._draw_time_series

    @property
    def draw_depth_profile(self) -> bool:
        return self._draw_depth_profile

    @staticmethod
    def read_plot_config(plot_name, plot_config, sources, variable_selections):
        if not isinstance(plot_config, dict):
            raise InvalidPlotConfig(
                'plot_config must be a dictionary that describe a plot.'
                'Received: {}'.format(str(plot_config))
            )

        for yaml_field in plot_config:
            if yaml_field not in PlotConfig.CONFIG_FIELDS:
                raise InvalidPlotConfig(
                    'Field "{}" is not an admissible field for a '
                    'plot_config'.format(yaml_field)
                )

        if 'source' not in plot_config:
            raise InvalidPlotConfig('Mandatory field "source" is missing')
        source_name = plot_config['source']

        if source_name not in sources:
            raise InvalidPlotConfig(
                'Plot "{}" uses source "{}", that has not been defined'.format(
                    plot_name,
                    source_name
                )
            )
        plot_source = sources[source_name]

        if 'variables' not in plot_config:
            raise InvalidPlotConfig('Mandatory field "variables" is missing')
        variable_selection_name = plot_config['variables']
        if variable_selection_name not in variable_selections:
            raise InvalidPlotConfig(
                'Plot "{}" uses "{}" variable selection, that has not been '
                'defined'.format(plot_name, variable_selection_name)
            )
        variables = variable_selections[variable_selection_name]

        if 'coast_index' not in plot_config:
            if plot_source.coast_index is not None:
                coast_index = plot_source.coast_index
            else:
                coast_index = DEFAULT_COAST_INDEX
        else:
            coast_index = int(plot_config['coast_index'])

        indicator = DEFAULT_INDICATOR
        if 'indicator' in plot_config:
            indicator = int(plot_config['indicator'])

        plot_builder = SingleLineInputData.builder(coast_index, indicator)

        plot_config_kwargs = {}

        def read_boolean(field_name):
            def bool_casting(x_val):
                if x_val is None:
                    x_val = False
                if x_val is not True and x_val is not False:
                    raise InvalidPlotConfig(
                        'The field "{}" must contain a boolean value: '
                        'received {}'.format(field_name, x_val)
                    )
                return x_val
            return bool_casting

        def add_to_kwargs(field_name, casting: Callable[[str], Any] = str):
            if field_name in plot_config:
                plot_config_kwargs[field_name] = casting(
                    plot_config[field_name]
                )

        add_to_kwargs('color')
        add_to_kwargs('zorder', int)
        add_to_kwargs('alpha', float)
        add_to_kwargs('alpha_for_time_series', float)
        add_to_kwargs('active', read_boolean('active'))
        add_to_kwargs('draw_depth_profile', read_boolean('draw_depth_profile'))
        add_to_kwargs('draw_time_series', read_boolean('draw_time_series'))

        # If the legend is not specified, is the name of the plot. If it is
        # "None" or "False", it is None
        if 'legend' in plot_config:
            if plot_config['legend'] is None:
                plot_config_kwargs['legend'] = None
            elif bool(plot_config['legend']) is False:
                plot_config_kwargs['legend'] = None
            else:
                plot_config_kwargs['legend'] = str(plot_config['legend'])
        else:
            plot_config_kwargs['legend'] = plot_name

        plot_filter = None
        if 'filter' in plot_config:
            f_description = plot_config['filter']
            if f_description is not None and f_description is not False:
                plot_filter = read_filter_description(f_description)

        return PlotConfig(
            plot_name,
            plot_source,
            variables,
            plot_builder,
            data_filter=plot_filter,
            **plot_config_kwargs
        )


@dataclass
class Config:
    plots: Tuple[PlotConfig, ...]
    time_series_options: TimeSeriesOptions
    depth_profiles_options: DepthProfilesOptions = DepthProfilesOptions()
    output_options: OutputOptions = OutputOptions()
    variable_labels: Mapping[str, str] = field(
        default_factory=lambda: defaultdict(dict)
    )


def read_config_from_file(config_path):
    with open(config_path, 'r') as f:
        current_config = read_config(f)
    return current_config


def read_config(config_datastream):
    yaml_content = yaml.safe_load(config_datastream)

    for yaml_field in yaml_content:
        if yaml_field not in EXPECTED_FIELDS:
            raise InvalidConfigFile(
                'Invalid field "{}" found inside config file. Admissible '
                'fields are: {}'.format(yaml_field, ', '.join(EXPECTED_FIELDS))
            )

    if 'sources' not in yaml_content:
        raise InvalidConfigFile(
            'Field "sources" not found inside the config file'
        )

    # Read the sources
    sources = {}
    sources_raw = yaml_content['sources']
    if not isinstance(sources_raw, dict):
        raise InvalidConfigFile(
            'Field "sources" must be a dictionary that associates a '
            'name to a path'
        )
    source_fields = {'path', 'meshmask', 'coast_index', 'mask_var_name'}
    for source_name, source_config in sources_raw.items():
        if not isinstance(source_config, dict):
            raise InvalidConfigFile(
                'A source must be associated to a dictionary that describes '
                'its path and its meshmask. Invalid field for source {}: '
                '{}'.format(source_name, source_config)
            )
        for yaml_field in source_config:
            if yaml_field not in source_fields:
                raise InvalidConfigFile(
                    'Invalid field "{}" for source "{}"'.format(
                        source_name,
                        yaml_field
                    )
                )
        if 'path' not in source_config:
            raise InvalidConfigFile(
                'No field "path" specified for source {}'.format(source_name)
            )
        if 'meshmask' not in source_config:
            raise InvalidConfigFile(
                'No field "meshmask" specified for source {}'.format(
                    source_name
                )
            )
        source_path = Path(str(source_config['path']))
        source_meshmask = Path(str(source_config['meshmask']))

        if 'coast_index' not in source_config:
            coast_index = None
        else:
            coast_index = int(source_config['coast_index'])

        if 'mask_var_name' in source_config:
            mask_var_name = str(source_config['mask_var_name'])
        else:
            mask_var_name = DEFAULT_MASK_VAR_NAME

        sources[str(source_name)] = DataDirSource(
            source_path,
            source_meshmask,
            coast_index,
            mask_var_name
        )

    if 'variable_selections' not in yaml_content:
        raise InvalidConfigFile(
            'Field "variable_selections" not found inside the config file'
        )

    # Read the selections of the variables
    variable_selections = {}
    variable_selections_raw = yaml_content['variable_selections']
    if not isinstance(variable_selections_raw, dict):
        raise InvalidConfigFile(
            'Field "variable_selections" must be a dictionary that associates '
            'a name to a list of biogeochemical variables'
        )
    for selection_name, selection_vars in variable_selections_raw.items():
        if not isinstance(selection_vars, Iterable):
            raise InvalidConfigFile(
                'Content of variable selection "{}" is not a valid '
                'list: {}'.format(selection_name, str(selection_vars))
            )
        variable_selections[str(selection_name)] = \
            tuple(str(i) for i in selection_vars)

    # Read the labels that we will use for variables
    if 'variable_labels' not in yaml_content:
        variable_labels = {}
    else:
        variable_labels_field = yaml_content['variable_labels']
        if not isinstance(variable_labels_field, Mapping):
            raise InvalidConfigFile(
                'The field variable_labels must be a dictionary that '
                'associates a variable name to its corresponding label'
            )
        variable_labels = {
            v: str(variable_labels_field[v]) for v in variable_labels_field
        }

    if 'time_series' in yaml_content:
        time_series_options = read_time_series_options(
            yaml_content['time_series']
        )
    else:
        raise InvalidConfigFile(
            'Field "time_series" not found inside the config file'
        )

    if 'depth_profiles' in yaml_content:
        depth_profiles_options = read_depth_profiles_options(
            yaml_content['depth_profiles']
        )
    else:
        depth_profiles_options = read_depth_profiles_options()

    # Read the output options
    if 'output' not in yaml_content:
        raise InvalidConfigFile(
            'Field "output" not found inside the config file'
        )
    output = yaml_content['output']
    if not isinstance(output, Mapping):
        raise InvalidConfigFile(
            'Field "output" must be a dictionary with several fields: '
            '{}...'.format(', '.join(sorted(list(OUTPUT_FIELDS))))
        )
    for key in output:
        if key not in OUTPUT_FIELDS:
            raise InvalidConfigFile(
                'Field "{}" is not allowed inside field "output"'.format(
                    key
                )
            )
    dpi = DEFAULT_DPI
    if 'dpi' in output:
        dpi = int(output['dpi'])
    output_dir = None
    if 'output_dir' in output:
        output_dir = read_output_dir(str(output['output_dir']))

    fig_size = DEFAULT_FIG_SIZE
    if 'fig_size' in output:
        try:
            fig_size = read_pair_of_integers(output['fig_size'])
        except ValueError:
            raise InvalidConfigFile('Invalid fig_size parameter')

    fig_ratio = DEFAULT_FIG_RATIO
    if 'fig_ratio' in output:
        try:
            fig_ratio_raw = read_pair_of_integers(output['fig_ratio'])
            fig_ratio_gcd = gcd(fig_ratio_raw[0], fig_ratio_raw[1])
            fig_ratio = (
                fig_ratio_raw[0] // fig_ratio_gcd,
                fig_ratio_raw[1] // fig_ratio_gcd
            )
        except ValueError:
            raise InvalidConfigFile('Invalid fig_ratio parameter')

    output_name = DEFAULT_OUTPUT_NAME
    if 'output_name' in output:
        output_name = str(output_name)

    output_options = OutputOptions(
        output_name=output_name,
        dpi=dpi,
        fig_size=fig_size,
        fig_ratio=fig_ratio,
        output_dir=output_dir
    )

    # Finally, the most complicated part; reading the plots
    if 'plots' not in yaml_content:
        raise InvalidConfigFile(
            'Field "plots" not found inside the config file'
        )
    plots = []
    plots_raw = yaml_content['plots']
    if not isinstance(plots_raw, dict):
        raise InvalidConfigFile(
            'Fields "plots" must be a dictionary that associates to a plot'
            'name its configuration'
        )
    for plot_name, plot_config in plots_raw.items():
        current_plot = PlotConfig.read_plot_config(
            plot_name,
            plot_config,
            sources,
            variable_selections
        )
        plots.append(current_plot)
    plots = tuple(p for p in plots if p.is_active())

    return Config(
        plots=plots,
        time_series_options=time_series_options,
        depth_profiles_options=depth_profiles_options,
        output_options=output_options,
        variable_labels=variable_labels
    )
