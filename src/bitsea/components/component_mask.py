from typing import Union

import numpy as np

from bitsea.commons.mask import Mask
from bitsea.components.find_component import find_component


class ComponentMask:
    """
    A class to identify and manage connected components in a mask.

    The mask can be a boolean array or a bit.sea `Mask` object. The class
    allows for decomposing the mask into its connected components, making it
    easy to retrieve individual components and analyze their properties.
    Components are assigned unique integer IDs in the output array.
    """

    def __init__(self, mask: Union[Mask, np.ndarray]):
        """
        Initializes the ComponentMask object by identifying connected
        components in the input mask.

        Args:
            mask: A boolean array or a bit.sea`Mask` object representing
                the input binary mask.
        """
        self._original_mask = mask

        work_area = self._original_mask.copy()
        n_components = 0

        self._components = np.full_like(work_area, -1, dtype=np.int32)
        self._n_cells_per_component = []

        while True:
            # Get the indices of the first element that is non-zero
            # (we exploit the fact that `True` is the biggest values of a
            # boolean array)
            first_true = np.unravel_index(np.argmax(work_area), work_area.shape)

            # If there are only `False` in the boolean, argmax should have
            # returned the (0, 0, ..., 0) because it is the first occurrence of
            # the maximum value (which is `False`). In this case, we are done
            if not work_area[first_true]:
                break

            # noinspection PyTypeChecker
            current_component = find_component(work_area, first_true)
            self._components[current_component] = n_components
            self._n_cells_per_component.append(
                np.count_nonzero(current_component)
            )
            work_area[current_component] = False
            n_components += 1

        self._n_components = n_components
        self._components.setflags(write=False)

    def get_component(self, component_id: int) -> np.ndarray:
        """
        Returns a boolean mask representing the specified connected component.

        Args:
            component_id: The ID of the connected component to retrieve. Should
                be between 0 and the total number of components - 1.

        Returns:
            A boolean array representing the mask of the specified component.

        Raises:
            IndexError: If the specified `component_id` is out of range.
        """
        if component_id < 0 or component_id >= self._n_components:
            raise IndexError(
                f"This mask has {self._n_components} components, numbered "
                f"between 0 and {self._n_components - 1}; unable to get "
                f"{component_id}"
            )
        return self._components == component_id

    @property
    def components(self):
        """
        Provides the full component mask with each connected component
        assigned a unique integer ID.

        The values that were `False` in the original mask are assigned -1. Each
        component is assigned a unique integer ID starting from 0.

        Returns:
            np.ndarray: An array of the same shape as the original mask, where
                each connected component is represented by its unique ID.
        """
        return self._components

    @property
    def n_components(self) -> int:
        """
        Returns the total number of connected components found in the input
        mask.

        Returns:
            int: The number of connected components in the mask.
        """
        return self._n_components

    def get_biggest_component(self) -> int:
        """
        Identifies the largest connected component in the input mask based on
        the number of cells.

        Returns:
            int: The ID of the largest connected component.
        """
        return int(np.argmax(self._n_cells_per_component))
