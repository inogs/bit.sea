# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import re

"""Helper functions"""

def is_number(val):
    """Tells if a value is a number type

    Args:
        - *val*: the value to test

    Returns:
        - True if val is a number type.
        - False if val is not a number type.
    """
    return isinstance(val, (int, long, float, complex, np.float32))

def get_date_string(s):
    """Finds a date in YYYYMMDD format into a bigger string.

    Args:
        - *s*: the string.

    Returns: a tuple with the date string in YYYY-MM-DD format as first element
    and YYYYMMDD format as second element. If no date has been found returns
    two empty strings.
    """
    m = re.search('.*([0-9]{4})([0-9]{2})([0-9]{2}).*',s)
    if m is None:
        return "", ""
    longdate = "%s-%s-%s" % (m.group(1), m.group(2), m.group(3))
    shortdate = "%s%s%s" % (m.group(1), m.group(2), m.group(3))
    return longdate, shortdate
