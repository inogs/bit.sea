from typing import overload
from typing import Union

import numpy as np


@overload
def is_leap(year: Union[int, np.int32, np.int64]):
    ...

def is_leap(year: np.ndarray) -> bool:
    """Check if a year is leap or not

    Args:
        year (int or array of ints): the year

    Returns:
         If year is an integer, returns a boolean value. Otherwise, returns
         an array of booleans, one for each year.
    """
    output = np.logical_and(
        year % 4 == 0,
        np.logical_or(year % 100 != 0, year % 400 == 0)
    )

    # If the result is not an array, cast it to a standard Python bool
    # (otherwise we get a numpy.bool_ object)
    if output.ndim == 0:
        output = bool(output)

    return output
