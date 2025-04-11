from itertools import product as cart_prod
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

    @staticmethod
    def _fill_holes(component: np.ndarray) -> np.ndarray:
        """
        Given a boolean mask, fills the holes that are within the component,
        returning a boolean mask that is True on all the points that are in
        the component or in all the points that are surrounded by the
        component.

        By definition, a "hole" is a region of `False` values that is not
        connected to the boundary of the mask.
        """
        # We allocate the space for the final output
        final_output = np.empty_like(component)
        final_output[:] = component

        # We build an object to check which ones of the components of `False`
        # are holes
        inverse_cm = ComponentMask(~component)
        for i in range(inverse_cm.n_components):
            current_component = inverse_cm.get_component(i)

            # Check if this component reaches the boundary; we need to build
            # all the slices of the form [:, :, ..., :, 0, :, ..., :] or
            # [:, :, ..., :, -1, :, ..., :] to analyze all the boundaries
            on_boundary = False
            for axis in range(component.ndim):
                for side in (0, -1):
                    boundary_slice = tuple(
                        side if k == axis else slice(None)
                        for k in range(component.ndim)
                    )
                    if np.count_nonzero(current_component[boundary_slice]) > 0:
                        on_boundary = True
                        break
                if on_boundary:
                    break

            # if we are not on boundary, this is a hole and must be added to
            # the output
            if not on_boundary:
                np.logical_or(final_output, current_component, out=final_output)

        return final_output

    @staticmethod
    def _get_boundary_cells(
        mask: np.ndarray, closed_boundary: bool = False
    ) -> np.ndarray:
        """Given a boolean mask, returns the cells that are inside the mask
        (i.e. that are `True`) and that are on the boundary of the mask.

        In other words, it returns the cells whose value is `True` but that
        have a neighbour that is `False`. With the word "neighbour", in this
        case, we also include cells that can be reached by going one step
        in every direction, including the diagonals.

        Args:
            mask: A boolean array representing the input mask.
            closed_boundary: If True, the cells inside the mask that are on the
                boundary of the domain will be also included in the output.
        """
        enlarged_data = np.pad(
            ~mask, 1, mode="constant", constant_values=closed_boundary
        )

        # To determine if a cell has a neighbour that is outside the mask, is
        # enough to shift the indices of tha mask by 1 on every possible
        # direction. If we apply the following slices to the enlarged_data we
        # may get the original values (slice associated to 0), the original
        # values moved towards the right (slice associated to 1) or ti the left
        # (slice associated to -1). By applying all the possible choice of the
        # slices for each axis, we get all the possible neighbours.
        slice_shifts = {0: slice(1, -1), 1: slice(2, None), -1: slice(None, -2)}

        # After the loop, this array will contain `True` if the point at that
        # position has a neighbour that is outside the mask
        neighbour_outside = np.zeros_like(mask)

        # For every axis, we decide where we need to shift
        for shift in cart_prod((-1, 0, 1), repeat=mask.ndim):
            # If all the shifts are 0 (do not move this axis), we get the
            # original values, and we can skip this check
            if tuple(set(shift)) == (0,):
                continue
            current_slice = tuple(slice_shifts[i] for i in shift)
            np.logical_or(
                neighbour_outside,
                enlarged_data[current_slice],
                out=neighbour_outside,
            )

        return np.logical_and(mask, neighbour_outside)

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

    def get_component_boundary(
        self, component_id: int, closed_boundary=False
    ) -> np.ndarray:
        """
        Returns the cells of the component with the specified ID that are on
        the boundary of the component, i.e. that have at least one neighbour
        that is not in the component.

        We consider as "neighbours" also points that can be reached by moving
        along the diagonals. This ensures that the corner of a concave mask
        is also included in the boundary.

        Args:
            component_id: The ID of the connected component to retrieve.
            closed_boundary: If True, the cells inside the mask that are on the
                boundary of the domain will be also included in the output.
        """
        component_mask = self.get_component(component_id)
        return self._get_boundary_cells(
            component_mask, closed_boundary=closed_boundary
        )

    def get_component_external_boundary(self, component_id: int) -> np.ndarray:
        """
        Returns the cells of the component with the specified ID that are on
        the external boundary of the component.

        This function is very close to `get_component_boundary`, but it does
        not return the cells that are on the boundary of a hole of the
        component.

        The algorithm is as follows: the holes of the component are filled and
        then the boundary cells are computed.

        Be aware that some points that are part of the external boundary may not
        be part of the original mask. This may happen if there is a hole that
        includes the vertex of a concave angle of the external boundary. In this
        case, when the hole is filled by the algorithm, the vertex becomes part
        of the boundary. If you don't want this behaviour, intersect the output
        of this function with the original mask.
        """
        component_mask = self.get_component(component_id)
        component_without_holes = self._fill_holes(component_mask)
        return self._get_boundary_cells(component_without_holes)
