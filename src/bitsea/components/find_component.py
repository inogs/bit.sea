from collections.abc import Iterable
from typing import Tuple
from typing import Union

import numpy as np
from numba import jit

from bitsea.commons.mask import Mask
from bitsea.components.neighbour_schemes import CrossNeighbours
from bitsea.components.neighbour_schemes import NeighbourScheme


def find_component(
    mask: Union[Mask, np.ndarray],
    starting_point: Iterable[int, ...],
    neighbour_scheme: NeighbourScheme = CrossNeighbours(),
):
    """
    Identify a connected component within a binary mask using a specified
    neighbor scheme.

    This function finds all points connected to a given starting point in a
    binary mask.
    Connectivity is determined by a specified neighbor scheme, which defines the
    relationships between adjacent points. The output is a boolean array where
    only the points in the connected component, starting from the given point,
    are `True`.

    This function can be used for applications such as detecting lakes or
    verifying whether a water mask contains separate basins. Additionally, it
    can identify small islands by providing the logical negation of a water
    mask as the input (where `True` values represent non-water cells).

    Arguments:
        mask: A boolean array (or Mask). Must have the same dimensions as the
            starting point.
        starting_point: The coordinates of the starting point in the mask.
            This point must correspond to a `True` value in the mask.
        neighbour_scheme: An object defining the neighborhood connectivity
            scheme (e.g., cross-neighbors or full neighbors). Defaults to
            CrossNeighbours().

    Raises:
        ValueError: If the starting point does not correspond to a `True` value
        in the mask.

    Returns:
        np.ndarray: A boolean array of the same shape as the input mask, where
            only the connected component starting from the given point is
            marked as `True`.
    """
    starting_point = tuple(t for t in starting_point)
    if len(starting_point) != mask.ndim:
        raise ValueError(
            f"Starting point {starting_point} must have the same"
            f"dimensionality as the mask (which has shape {mask.shape})"
        )

    if not mask[starting_point]:
        raise ValueError('Starting point is not set to "True"')

    # `neighbours` is a 2D array where each row `i` represents a vector.
    # Adding this vector to the starting position gives the coordinates
    # of a neighboring point.
    neighbours = neighbour_scheme.get_indices(mask.ndim).astype(np.int32)

    # We keep the mask as a pure array, so it can be used by Numba
    mask_array = np.asarray(mask)

    # And here we delegate the function to Numba
    component_mask = _find_component_internal(
        mask_array.ravel(order="C"),
        np.asarray(mask.shape, dtype=np.int32),
        starting_point,
        neighbours,
    )
    return component_mask.astype(bool).reshape(mask.shape, order="C")


@jit(nopython=True)
def _find_component_internal(
    mask: np.ndarray,
    axes_limits: np.ndarray,
    starting_point: Tuple[int, ...],
    neighbours: np.ndarray,
):
    """
    Internal implementation of the find_component function that uses Numba for
    performance optimization.

    This function performs the connected component search based on the provided
    mask, starting point, and neighborhood scheme. It is optimized with
    Numba's JIT compilation.
    """
    dims = axes_limits.size
    n_neighbours = neighbours.shape[0]

    starting_point = np.array(starting_point, dtype=np.int32)
    # If the starting point contains negative indices, convert them to their
    # equivalent positive indices within the array dimensions.
    for k in range(dims):
        if starting_point[k] < 0:
            starting_point[k] += axes_limits[k]

    # Each element in the mask is identified by a single index. The
    # `index_transformation` array is used to map multidimensional
    # coordinates to their corresponding flat (1D) index in the array,
    # based on lexicographical order (starting from the last dimension).
    # The scalar product between the index_transformation array and an array
    # with the multidimensional coordinates returns the single index.
    index_transformation = np.ones(dims, dtype=np.int32)
    for i in range(1, dims):
        index_transformation[-i - 1] = (
            index_transformation[-i] * axes_limits[-i]
        )

    # An array to store the coordinates of the currently visited point.
    # This will be reused in the main loop for efficiency.
    current_coords = np.empty((dims,), dtype=np.int32)

    # This array stores points that have been found as part of the component
    # and need to be visited to examine their neighbors. Since Numba does
    # not support dynamic lists, a large array is pre-allocated and initialized
    # with -1 (indicating "empty" slots).
    to_be_visited = np.full(mask.size, -1, dtype=np.int32)

    # Mark the starting point as the first point to visit.
    # Since Numba lacks a `.dot` method for integers, the index
    # of the starting point is computed manually using a loop.
    to_be_visited[0] = 0
    for i in range(dims):
        to_be_visited[0] += starting_point[i] * index_transformation[i]

    # This is the final output. It is True on the points that can be reached
    # starting from the starting_point and moving following the neighbours
    # arrays
    component_mask = np.zeros((mask.size,), dtype=np.bool_)
    component_mask[to_be_visited[0]] = True

    # How many points to we have still to visit (1, the starting point)
    visiting_index = 1

    while visiting_index > 0:
        visiting_index -= 1
        currently_visited = to_be_visited[visiting_index]
        assert currently_visited >= 0, (
            "We should not be visiting an empty point"
        )

        # Fill the current_coords vector with the coordinates of the point we
        # are visiting
        current_coords[0] = currently_visited // index_transformation[0]
        for k in range(1, dims):
            current_coords[k] = (
                currently_visited % index_transformation[k - 1]
            ) // index_transformation[k]

        # Now we check what happens to each of its neighbours
        for n in range(n_neighbours):
            neighbour = neighbours[n]
            neighbour_coords = current_coords + neighbour

            # Maybe, if we move following the direction described by the
            # neighbour vector, we end up outside our array (for example,
            # because on the right boundary a cell does not have a neighbour
            # on its right). In this case, we ignore this point
            outside_mask = False
            for i in range(dims):
                if (
                    neighbour_coords[i] < 0
                    or neighbour_coords[i] >= axes_limits[i]
                ):
                    outside_mask = True
                    break
            if outside_mask:
                continue

            # Here we compute the index of the neighbour cell
            neighbour_index = 0
            for i in range(dims):
                neighbour_index += neighbour_coords[i] * index_transformation[i]

            # Point on the land or outside the mask
            if not mask[neighbour_index]:
                continue

            # Point already visited
            if component_mask[neighbour_index]:
                continue

            # Ok, we have found a new element for this component; we add it to
            # our mask
            component_mask[neighbour_index] = True

            # Let us remember that we need to visit this new point to check the
            # status of its neighbours
            to_be_visited[visiting_index] = neighbour_index
            visiting_index += 1

    return component_mask
