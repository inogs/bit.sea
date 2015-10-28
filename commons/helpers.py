# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np

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

