from dataclasses import dataclass, field, InitVar
from enum import Enum
from types import MappingProxyType
from typing import Any, Dict, Mapping, Tuple, Union
import re


class InvalidAlgorithmDescription(Exception):
    pass


class DepthProfileAlgorithm(Enum):
    STANDARD = 1
    SEASONAL = 2
    EVOLUTION = 3


_ALGORITHM_ASSOCIATION = {
    'standard': DepthProfileAlgorithm.STANDARD,
    'seasonal': DepthProfileAlgorithm.SEASONAL,
    'evolution': DepthProfileAlgorithm.EVOLUTION
}


@dataclass(frozen=True)
class DepthProfileMode:
    algorithm: DepthProfileAlgorithm = DepthProfileAlgorithm.STANDARD
    config_: InitVar = None
    config: Mapping[str, Any] = field(init=False)

    def __post_init__(self, config_: Union[Mapping[str, Any], None] = None):
        if config_ is None:
            _config = {}
        super().__setattr__(
            '_config',
            {a: b for a, b in config_.items()}
        )
        super().__setattr__('config', MappingProxyType(config_))


DEFAULT_DEPTH_PROFILE_MODE = DepthProfileMode(
    algorithm=DepthProfileAlgorithm.STANDARD,
    config_={}
)


MODE_DESCRIPTION_MASK = re.compile(
    r'^\s*(?P<algorithm_name>[a-zA-Z][a-zA-Z0-9_]*)((?P<args>\(.*\))?)\s*$'
)


def read_depth_profile_mode(read_profile_raw: str) -> DepthProfileMode:
    mask_match = MODE_DESCRIPTION_MASK.match(read_profile_raw)
    if mask_match is None:
        raise InvalidAlgorithmDescription(
            'Invalid algorithm description! An algorithm must satisfy the '
            'following regular expression: {}\nReceived: "{}"'.format(
                MODE_DESCRIPTION_MASK,
                read_profile_raw
            )
        )
    algorithm_raw = mask_match.group('algorithm_name')
    raw_options = mask_match.group('args')

    if algorithm_raw.lower() not in _ALGORITHM_ASSOCIATION:
        raise InvalidAlgorithmDescription(
            'Unknown algorithm "{}"; the only accepted values are: {}'.format(
                algorithm_raw,
                ','.join(['"' + k + '"' for k in _ALGORITHM_ASSOCIATION.keys()])
            )
        )
    algorithm = _ALGORITHM_ASSOCIATION[algorithm_raw.lower()]

    if raw_options is not None:
        # Remove the parenthesis
        raw_options = raw_options[1:-1].strip()

        if raw_options == '':
            raw_options = tuple()
        else:
            raw_options = tuple(k.strip() for k in raw_options.split(','))
    else:
        raw_options = tuple()

    options: Dict[str, Any] = {}

    if algorithm is DepthProfileAlgorithm.STANDARD and len(raw_options) != 0:
        raise InvalidAlgorithmDescription(
            'Standard algorithm does not accept arguments: received '
            '"{}"'.format(raw_options)
        )

    if algorithm is DepthProfileAlgorithm.SEASONAL and len(raw_options) != 0:
        if len(raw_options) > 2:
            raise InvalidAlgorithmDescription(
                'Seasonal algorithm accepts only one argument, the "mode"; '
                'received {} arguments: {}'.format(
                    len(raw_options),
                    raw_options
                )
            )
        option = raw_options[0]
        if '=' in option:
            key = option.split('=')[0].strip()
            option = '='.join(option.split('=')[1:])
            if key.lower() != 'mode':
                raise InvalidAlgorithmDescription(
                    'The only argument accepted by a seasonal algorithm is '
                    'the mode; received "{}"'.format(key)
                )
        option = option.strip()
        if option.lower() not in ('square', 'inline'):
            raise InvalidAlgorithmDescription(
                'The only accepted arguments for a seasonal algorithm are '
                '"square", "inline"; received {}'.format(option)
            )
        options = {'mode': option}

    if algorithm is DepthProfileAlgorithm.SEASONAL and len(raw_options) == 0:
        options = {'mode': 'square'}

    if algorithm is DepthProfileAlgorithm.EVOLUTION:
        if len(raw_options) == 0:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires an integer as an argument; '
                'none received'
            )
        if len(raw_options) > 1:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires a single integer as an '
                'argument; received {} arguments: {}'.format(
                    len(raw_options),
                    raw_options
                )
            )

        try:
            i_value = int(raw_options[0])
        except ValueError:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires an integer as an argument: '
                'received "{}"'.format(raw_options[0])
            )

        if i_value < 2:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires an integer number greater than'
                '1: received {}'.format(i_value)
            )

        options = {'plots': i_value}

    return DepthProfileMode(algorithm, options)


def get_depth_profile_plot_grid(depth_profile_mode: DepthProfileMode) \
        -> Tuple[int, int]:
    algorithm = depth_profile_mode.algorithm
    if algorithm == DepthProfileAlgorithm.STANDARD:
        return 1, 1
    elif algorithm == DepthProfileAlgorithm.EVOLUTION:
        return 1, 1
    elif algorithm == DepthProfileAlgorithm.SEASONAL:
        if depth_profile_mode.config['mode'] == 'square':
            return 2, 2
        elif depth_profile_mode.config['mode'] == 'inline':
            return 1, 4
        else:
            raise ValueError(
                'Invalid mode "{}" for a seasonal algorithm'.format(
                    depth_profile_mode.config['mode']
                )
            )
    raise ValueError('Unknown algorithm: {}'.format(algorithm))
