from collections import namedtuple
from enum import Enum
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


DepthProfileMode = namedtuple(
    'DepthProfileMode',
    ('algorithm', 'config')
)


DEFAULT_DEPTH_PROFILE_MODE = (DepthProfileAlgorithm.STANDARD, ())


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
        raise InvalidAlgorithmDescription(
            'Seasonal algorithm does not accept arguments: received '
            '"{}"'.format(options)
        )

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
