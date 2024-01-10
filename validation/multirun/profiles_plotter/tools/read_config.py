from collections import namedtuple
from collections.abc import Iterable
import yaml
from pathlib import Path
from sys import version_info
from typing import Any

from tools.data_object import PickleDataObject
from tools.depth_profile_algorithms import DEFAULT_DEPTH_PROFILE_MODE, \
    read_depth_profile_mode
from filters.read_filter_description import read_filter_description

if version_info[1] < 9:
    from typing import Callable
else:
    from collections.abc import Callable


EXPECTED_FIELDS = {
    'sources',
    'variable_selections',
    'plots',
    'levels',
    'depth_profile_mode',
    'output_dir'
}


class InvalidConfigFile(Exception):
    pass


class InvalidPlotConfig(InvalidConfigFile):
    pass


Config = namedtuple(
    'Config',
    ('plots', 'levels', 'depth_profile_mode', 'output_dir')
)

Source = namedtuple('Source', ('path', 'meshmask'))


def read_number(number_str):
    try:
        return int(number_str)
    except ValueError:
        return float(number_str)


class PlotConfig:
    CONFIG_FIELDS = {
        'source',
        'variables',
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

    def __init__(self, name, source, variables, data_filter=None, active=True,
                 color=None, alpha=None, alpha_for_time_series=None,
                 legend=None, zorder=None, draw_time_series=True,
                 draw_depth_profile=True):
        if name in self.PLOT_NAMES:
            raise ValueError(
                'A plot named "{}" has already been created'.format(name)
            )
        self.PLOT_NAMES.add(name)

        self.name = name
        self._source = source
        self.variables = variables

        self._active = active
        self._filter = data_filter

        self._color = color
        self._alpha = alpha
        self._alpha_for_time_series = alpha_for_time_series
        self._legend = legend
        self._zorder = zorder

        self._draw_time_series = bool(draw_time_series)
        self._draw_depth_profile = bool(draw_depth_profile)

    def get_plot_data(self, variable):
        if variable not in self.variables:
            raise ValueError(
                'Plot "{}" does not contain variable "{}"'.format(
                    self.name,
                    variable
                )
            )

        data_object = PickleDataObject(self.source.path, variable)
        if self._filter is None:
            return data_object

        return self._filter.get_filtered_object(data_object)

    @property
    def source(self):
        return self._source

    def is_active(self):
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
    def legend(self):
        return self.name if self._legend is None else self._legend

    @property
    def zorder(self):
        return self._zorder

    @property
    def draw_time_series(self):
        return self._draw_time_series

    @property
    def draw_depth_profile(self):
        return self._draw_depth_profile

    @staticmethod
    def read_plot_config(plot_name, plot_config, sources, variable_selections):
        if not isinstance(plot_config, dict):
            raise InvalidPlotConfig(
                'plot_config must be a dictionary that describe a plot.'
                'Received: {}'.format(str(plot_config))
            )

        for field in plot_config:
            if field not in PlotConfig.CONFIG_FIELDS:
                raise InvalidPlotConfig(
                    'Field "{}" is not an admissible field for a '
                    'plot_config'.format(field)
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
        add_to_kwargs('legend')
        add_to_kwargs('zorder', int)
        add_to_kwargs('alpha', float)
        add_to_kwargs('alpha_for_time_series', float)
        add_to_kwargs('active', read_boolean('active'))
        add_to_kwargs('draw_depth_profile', read_boolean('draw_depth_profile'))
        add_to_kwargs('draw_time_series', read_boolean('draw_time_series'))

        plot_filter = None
        if 'filter' in plot_config:
            f_description = plot_config['filter']
            if f_description is not None and f_description is not False:
                plot_filter = read_filter_description(f_description)

        return PlotConfig(
            plot_name,
            plot_source,
            variables,
            data_filter=plot_filter,
            **plot_config_kwargs
        )


def read_config_from_file(config_path):
    with open(config_path, 'r') as f:
        current_config = read_config(f)
    return current_config


def read_config(config_datastream):
    yaml_content = yaml.safe_load(config_datastream)

    for field in yaml_content:
        if field not in EXPECTED_FIELDS:
            raise InvalidConfigFile(
                'Invalid field "{}" found inside config file. Admissible '
                'fields are: {}'.format(field, ', '.join(EXPECTED_FIELDS))
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
    for source_name, source_config in sources_raw.items():
        if not isinstance(source_config, dict):
            raise InvalidConfigFile(
                'A source must be associated to a dictionary that describes '
                'its path and its meshmask. Invalid field for source {}: '
                '{}'.format(source_name, source_config)
            )
        for field in source_config:
            if field not in ('path', 'meshmask'):
                raise InvalidConfigFile(
                    'Invalid field "{}" for source "{}"'.format(
                        source_name,
                        field
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
        sources[str(source_name)] = Source(source_path, source_meshmask)

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

    # Read levels, ensuring that they are a list of integers or floats
    if 'levels' not in yaml_content:
        raise InvalidConfigFile(
            'Field "levels" not found inside the config file'
        )
    levels_raw = yaml_content['levels']
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

    if 'depth_profile_mode' in yaml_content:
        depth_profile_mode = read_depth_profile_mode(
            yaml_content['depth_profile_mode']
        )
    else:
        depth_profile_mode = DEFAULT_DEPTH_PROFILE_MODE

    # Read the output dir
    if 'output_dir' not in yaml_content:
        raise InvalidConfigFile(
            'Field "output_dir" not found inside the config file'
        )
    output_dir = Path(str(yaml_content['output_dir']))

    # Finally, the most complicated part; reading the plots
    if 'plots' not in yaml_content:
        raise InvalidConfigFile(
            'Field "plots" not found inside the config file'
        )
    plots = []
    plots_raw = yaml_content['plots']
    if not isinstance(plots_raw, dict):
        raise InvalidConfigFile(
            'Fields "plots" must be a dictionary tat associates to a plot name'
            'its configuration'
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
        levels=levels,
        depth_profile_mode=depth_profile_mode,
        output_dir=output_dir
    )
