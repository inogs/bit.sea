from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Tuple
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


@dataclass
class DepthProfileMode:
    algorithm: DepthProfileAlgorithm = DepthProfileAlgorithm.STANDARD
    config: Mapping[str, Any] = field(
        default_factory=lambda: defaultdict(dict)
    )


DEFAULT_DEPTH_PROFILE_MODE = DepthProfileMode(
    algorithm=DepthProfileAlgorithm.STANDARD,
    config={}
)


MODE_DESCRIPTION_MASK = re.compile(
    r'^\s*(?P<algorithm_name>[a-zA-Z][a-zA-Z0-9_]*)((?P<args>\(.*\))?)\s*$'
)


def read_depth_profile_mode(read_profile_raw):
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
    options = mask_match.group('args')

    if algorithm_raw.lower() not in _ALGORITHM_ASSOCIATION:
        raise InvalidAlgorithmDescription(
            'Unknown algorithm "{}"; the only accepted values are: {}'.format(
                algorithm_raw,
                ','.join(['"' + k + '"' for k in _ALGORITHM_ASSOCIATION.keys()])
            )
        )
    algorithm = _ALGORITHM_ASSOCIATION[algorithm_raw.lower()]

    if options is not None:
        # Remove the parenthesis
        options = options[1:-1].strip()

        if options == '':
            options = tuple()
        else:
            options = tuple(k.strip() for k in options.split(','))
    else:
        options = tuple()

    if algorithm is DepthProfileAlgorithm.STANDARD and len(options) != 0:
        raise InvalidAlgorithmDescription(
            'Standard algorithm does not accept arguments: received '
            '"{}"'.format(options)
        )

    if algorithm is DepthProfileAlgorithm.SEASONAL and len(options) != 0:
        if len(options) > 2:
            raise InvalidAlgorithmDescription(
                'Seasonal algorithm accepts only one argument, the "mode"; '
                'received {} arguments: {}'.format(len(options), options)
            )
        option = options[0]
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

    if algorithm is DepthProfileAlgorithm.SEASONAL and len(options) == 0:
        options = {'mode': 'square'}

    if algorithm is DepthProfileAlgorithm.EVOLUTION:
        if len(options) == 0:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires an integer as an argument; '
                'none received'
            )
        if len(options) > 1:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires a single integer as an '
                'argument; received {} arguments: {}'.format(
                    len(options),
                    options
                )
            )

        try:
            i_value = int(options[0])
        except ValueError:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires an integer as an argument: '
                'received "{}"'.format(options[0])
            )

        if i_value < 2:
            raise InvalidAlgorithmDescription(
                '"Evolution" algorithm requires an integer number greater than'
                '1: received {}'.format(i_value)
            )

        options = (i_value,)

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
