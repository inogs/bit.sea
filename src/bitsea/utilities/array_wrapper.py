from typing import Union

import numpy as np
from numpy.typing import DTypeLike


class ArrayWrapper:
    """
    `ArrayWrapper` is a protocol that can be inherited by other classes if they
    need to contain a NumPy array and expose most of its methods to the user.

    The goal is for the resulting object to behave like a NumPy array in most
    functions. This is achieved by implementing the `__array__` protocol from
    NumPy.

    To access the underlying wrapped array directly, the `as_array` and
    `as_mutable_array` methods are provided.

    When retrieving data from this object using the `getitem` method (i.e.,
    square brackets), the returned slice will have the `write` flag set to
    `False`, ensuring the data of this object remains unmodified.
    """

    def __init__(self, *, wrapped_data: np.ndarray, **kwargs):
        self._data_array = np.asarray(wrapped_data).view()
        super().__init__(**kwargs)

    def as_array(self) -> np.ndarray:
        """
        Retrieves the NumPy array wrapped by this object and returns it.
        The returned array is a view of the original data, with the writable
        flag set to `False` to prevent modifications to the data of this object.

        Returns:
            A NumPy array view of the original data with write access disabled.
        """
        output = self._data_array.view()
        output.setflags(write=False)
        return output

    def as_mutable_array(self):
        """
        Retrieves the NumPy array wrapped by this object and returns it.
        Unlike `as_array()`, the returned view is writable,
        allowing modifications to the data. If the returned array is
        modified, also the content of this object changes accordingly.

        Returns:
            A NumPy array view of the original data.
        """
        return self._data_array.view()

    def __array__(
        self,
        dtype: Union[DTypeLike, None] = None,
        copy: Union[bool, None] = None,
    ) -> np.ndarray:
        dtype = np.dtype(dtype) if dtype is not None else None

        copy_needed = False
        if dtype is not None:
            if self._data_array.dtype != dtype:
                copy_needed = True

        if copy_needed and copy is False:
            raise ValueError(
                "Unable to avoid copy while creating an array as requested.\n"
                "If using `np.array(obj, copy=False)` replace it with "
                "`np.asarray(obj)` to allow a copy when needed (no behavior "
                "change in NumPy 1.x).\n"
                "For more details, see "
                "https://numpy.org/devdocs/numpy_2_0_migration_guide.html"
                "#adapting-to-changes-in-the-copy-keyword."
            )

        if copy:
            copy_needed = True

        if not copy_needed:
            return self.as_array()

        return np.array(self._data_array, dtype=dtype, copy=True)

    @property
    def dtype(self):
        return self._data_array.dtype

    @property
    def ndim(self) -> int:
        return self._data_array.ndim

    @property
    def size(self):
        return self._data_array.size

    @property
    def shape(self):
        return self._data_array.shape

    def __getitem__(self, item) -> np.ndarray:
        output = self._data_array.__getitem__(item)

        # If `item` triggers a copy (i.e., if it uses advanced indices),
        # then we can return the selection as is. Otherwise, we set
        # the `write` flag to `False` to preserve the data of this
        # object
        if np.may_share_memory(output, self._data_array):
            output.setflags(write=False)

        return output

    def __setitem__(self, key, value):
        return self._data_array.__setitem__(key, value)


class BooleanArrayWrapper(ArrayWrapper):
    """
    A specialized version of `ArrayWrapper` that implements methods specific
    to NumPy arrays of booleans.
    """

    def __init__(self, *, wrapped_data: np.ndarray, **kwargs):
        if np.asarray(wrapped_data).dtype not in (bool, np.bool_):
            raise ValueError("Data must be a boolean array")

        super().__init__(wrapped_data=wrapped_data, **kwargs)

    def __invert__(self) -> np.ndarray:
        return np.logical_not(self._data_array)
