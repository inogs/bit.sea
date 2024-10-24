import hypothesis.strategies as st
import numpy as np
import pytest
from hypothesis import assume
from hypothesis import given

from bitsea.utilities.array_wrapper import ArrayWrapper
from bitsea.utilities.array_wrapper import BooleanArrayWrapper


@given(shape=st.lists(st.integers(min_value=1, max_value=10), max_size=4))
def test_array_wrapper_init(shape):
    shape = tuple(shape)
    test_cases = [
        np.random.choice([True, False], size=shape),
        np.random.random(size=shape),
        np.random.randint(low=-1000, high=1000, size=shape, dtype=int),
        np.random.randint(low=-1000, high=1000, size=shape, dtype=np.int32),
    ]

    for test_case in test_cases:
        array_wrapper = ArrayWrapper(test_case)

        assert test_case.shape == array_wrapper.shape
        assert test_case.dtype == array_wrapper.dtype

        assert np.allclose(test_case, array_wrapper)


@given(shape=st.lists(st.integers(min_value=1, max_value=10), max_size=4))
def test_array_wrapper_ndim_method(shape):
    shape = tuple(shape)

    test_cases = [
        np.random.choice([True, False], size=shape),
        np.random.random(size=shape),
        np.random.randint(low=-1000, high=1000, size=shape, dtype=int),
        np.random.randint(low=-1000, high=1000, size=shape, dtype=np.int32),
    ]

    for test_case in test_cases:
        array_wrapper = ArrayWrapper(test_case)
        assert array_wrapper.ndim == len(shape)


@given(slice_i=st.slices(10), slice_j=st.slices(11))
def test_array_wrapper_get_item(slice_i, slice_j):
    test_data = np.linspace(0, 1, 10) * np.linspace(0, 1, 11)[:, np.newaxis]
    array_wrapper = ArrayWrapper(test_data)

    assert np.allclose(
        array_wrapper[slice_i, slice_j], test_data[slice_i, slice_j]
    )


@given(test_slice=st.slices(100))
def test_array_wrapper_get_item_is_read_only(test_slice):
    array_wrapper = ArrayWrapper(np.random.random(100))

    test_data = array_wrapper[test_slice]

    assume(test_data.size > 0)

    with pytest.raises(ValueError):
        test_data[:] = 0.0


@given(test_slice=st.slices(100))
def test_array_wrapper_set_item(test_slice):
    array_wrapper = ArrayWrapper(np.random.random(100))

    array_wrapper[test_slice] = 7
    assert np.allclose(array_wrapper[test_slice], 7)


@given(test_slice=st.slices(100))
def test_array_wrapper_as_mutable_array(test_slice):
    array_wrapper = ArrayWrapper(np.random.random(100))
    wrapper_data = array_wrapper.as_mutable_array()

    wrapper_data[test_slice] = 7

    assert np.allclose(array_wrapper[test_slice], 7)


@given(shape=st.lists(st.integers(min_value=1, max_value=10), max_size=4))
def test_boolean_array_wrapper_requires_boolean_arrays(shape):
    wrong_argument = np.random.random(size=shape)

    good_argument = np.random.choice([True, False], size=shape)

    with pytest.raises(ValueError):
        BooleanArrayWrapper(wrong_argument)

    assert BooleanArrayWrapper(good_argument).dtype == bool


@given(shape=st.lists(st.integers(min_value=1, max_value=10), max_size=4))
def test_boolean_array_wrapper_as_mask(shape):
    mask = BooleanArrayWrapper(np.random.choice([True, False], size=shape))

    test_data = np.full(fill_value=3, shape=shape, dtype=int)

    assert np.sum(test_data[mask]) == np.count_nonzero(mask) * 3
    assert np.sum(test_data[~mask]) == (mask.size - np.count_nonzero(mask)) * 3
