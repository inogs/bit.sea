import re
from filters.filter import ComposedFilter, InvalidFilterDescription
from filters.crop_time import TimeCroppingFilter
from filters.moving_average_filter import MovingAverageFilter


FILTER_DESCRIPTION_MASK = re.compile(
    r'^\s*(?P<filter_name>[a-zA-Z][a-zA-Z0-9_]*)((?P<args>\(.*\))?)\s*$'
)


# Add here the name of the filter with its association every time you introduce
# a new filter
_FILTER_ASSOCIATION = {
    'movingaverage': MovingAverageFilter,
    'croptime': TimeCroppingFilter
}


def read_filter_description(filter_description: str):
    if ';' in filter_description:
        # Multiple filters have been submitted: we analyze them one by one
        filters = tuple(
            read_filter_description(k) for k in filter_description.split(';')
        )
        return ComposedFilter(filters)

    # Now we analyze the description of a single filter
    mask_match = FILTER_DESCRIPTION_MASK.match(filter_description)
    if mask_match is None:
        raise InvalidFilterDescription(
            'Invalid filter description! A filter must satisfy the following '
            'regular expression: {}\nReceived: "{}"'.format(
                FILTER_DESCRIPTION_MASK,
                filter_description
            )
        )
    filter_name = mask_match.group('filter_name')
    filter_args = mask_match.group('args')

    filter_name_lower = filter_name.lower()
    if filter_name_lower in _FILTER_ASSOCIATION:
        # Remove the parenthesis
        if filter_args is not None:
            filter_args = filter_args[1:-1]
        return _FILTER_ASSOCIATION[filter_name_lower].initialize_from_string(
            filter_args
        )

    raise InvalidFilterDescription(
        'Invalid filter specified: "{}"\nAllowed values are: {}'.format(
            filter_name,
            ', '.join(sorted(list(_FILTER_ASSOCIATION.keys())))
        )
    )
