import numpy as np

from bitsea.utilities.time import is_leap


def test_is_leap():
    assert is_leap(2020) is True
    assert is_leap(2021) is False
    assert is_leap(1900) is False
    assert is_leap(2000) is True


def test_is_leap_vectorize():
    a = np.array([2020, 2021, 1900, 2000], dtype=np.int32)
    b = np.array([1999, 2020, 2021, 1900, 2000], dtype=int)

    a_leap = is_leap(a)
    b_leap = is_leap(b)

    for i, year in enumerate(a):
        assert a_leap[i] == is_leap(year)

    for i, year in enumerate(b):
        assert b_leap[i] == is_leap(year)
