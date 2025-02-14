"""
This module defines the framework and implementations for neighborhood schemes
used to find neighboring indices in multidimensional grids. It provides
abstract and concrete classes for neighborhood structures such as cross
neighbors and square-around neighbors.

Classes:
    NeighbourScheme: An abstract base class for defining neighborhood schemes.
    CrossNeighbours: Implements a neighborhood scheme including only orthogonal
        neighbors.
    SquareAroundNeighbours: Implements a neighborhood scheme including all
        neighbors in a square/cube.
"""

from abc import ABC
from abc import abstractmethod
from itertools import product as cart_product

import numpy as np


class NeighbourScheme(ABC):
    """
    An abstract base class for defining neighborhood schemes in a grid.

    A neighborhood scheme defines the relative positions (indices) of neighbors
    for a given central point in a multidimensional grid. Subclasses must
    provide concrete implementations of the `get_indices` method.
    """

    @abstractmethod
    def get_indices(self, dimension: int = 2) -> np.ndarray:
        """
        Get the relative indices for neighbors in the defined scheme.

        Parameters:
            dimension (int): The dimensionality of the grid (default is 2).

        Returns:
            np.ndarray: A 2D array where each row represents the relative
            coordinates of a neighboring point.
        """
        raise NotImplementedError


class CrossNeighbours(NeighbourScheme):
    """
    A neighborhood scheme that includes only orthogonal neighbors.

    In this scheme, at most two neighbors exist for each axis: one in the
    positive direction and one in the negative direction.
    """

    def get_indices(self, dimension: int = 2) -> np.ndarray:
        """
        Get the relative indices for cross neighbors in the specified dimension.

        Parameters:
            dimension (int): The dimensionality of the grid (default is 2).

        Returns:
            np.ndarray: A 2D array where each row represents the relative
            coordinates of a cross neighbor.
        """
        indices = tuple([] for _ in range(dimension))
        for d in range(dimension):
            for k in (-1, 1):
                for j in range(dimension):
                    if j == d:
                        indices[j].append(k)
                    else:
                        indices[j].append(0)
        return np.array(list(zip(*indices)))


class SquareAroundNeighbours(NeighbourScheme):
    """
    A neighborhood scheme that includes all neighbors in a square or cube.

    In this scheme, all neighbors in a given range (inclusive of diagonals)
    are included, except the center point itself.
    """

    def get_indices(self, dimension: int = 2):
        """
        Get the relative indices for square-around neighbors in the specified
        dimension.

        Parameters:
            dimension (int): The dimensionality of the grid (default is 2).

        Returns:
            np.ndarray: A 2D array where each row represents the relative
            coordinates of a square-around neighbor.
        """
        indices = tuple([] for _ in range(dimension))
        zero_coords = tuple(0 for _ in range(dimension))
        choices = (-1, 0, 1)
        for d in cart_product(choices, repeat=dimension):
            if d == zero_coords:
                continue
            for i in range(dimension):
                indices[i].append(d[i])
        return np.array(list(zip(*indices)))
