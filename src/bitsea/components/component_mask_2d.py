import sys
from itertools import product as cart_prod
from typing import Dict
from typing import List
from typing import Set
from typing import Tuple

import numpy as np

from bitsea.components.component_mask import ComponentMask


# All this code is just to define the class Point as a tuple of two integers;
# we need the "if" because TypeAlias exists only in Python 3.9
if sys.version_info >= (3, 10):
    from typing import TypeAlias

    Point: TypeAlias = Tuple[int, int]
else:
    Point = Tuple[int, int]


class NonEulerianPathError(ValueError):
    pass


class _PointGraph:
    """A graph representation of points and their connections on a 2D grid.

    This internal class manages points and their relationships in a 2D grid
    system.
    Points are represented as (x, y) integer coordinates and can be connected
    to other points via edges. The primary purpose is to model faces on a 2D
    grid, where points serve as vertices and edges represent the faces'
    boundaries.

    Note:
        This is an internal implementation class and should not be used
        directly by external users.

    Attributes:
        _points (Dict[Point, Set[Point]]): Dictionary storing points and their
            connected neighbors. Each point maps to a set of points it's
            connected to.
    """

    def __init__(self):
        self._points: Dict[Point, Set[Point]] = {}

    def add_point(self, point: Point) -> None:
        """Adds a point to the graph.

        If the point already exists, it does nothing."""
        if point not in self._points:
            self._points[point] = set()

    def add_edge(self, point1: Point, point2: Point) -> None:
        """Add an edge between two points.

        If the two points are not already part of the graph, they are added.
        If the edge already exists, it does nothing.
        """
        self.add_point(point1)
        self.add_point(point2)
        self._points[point1].add(point2)
        self._points[point2].add(point1)

    def __len__(self):
        return len(self._points)

    @property
    def n_points(self):
        return len(self)

    @property
    def n_edges(self):
        return sum(len(p) for p in self._points.values()) // 2

    def _get_connected_components(self) -> List[Set[Point]]:
        """Splits the graph points into disjoint connected components.

        A connected component is a set of points where any two points are
        connected by a path of edges. Points in different components have no
        path between them.

        Returns:
            List[Set[Point]]: A list of sets where:
                - Each set contains points that form a connected component
                - The union of all sets equals the graph's points
                - The intersection of any two sets is empty
                - Points in different sets have no connecting edges
                - Points in the same set are connected by a path of edges
        """
        connected_components: List[Set[Point]] = []
        missing = set(self._points.keys())

        while len(missing) > 0:
            current_component = set()
            to_visit = [missing.pop()]
            while len(to_visit) > 0:
                current_point = to_visit.pop()
                current_component.add(current_point)
                to_visit.extend(self._points[current_point] - current_component)
            missing -= current_component
            connected_components.append(current_component)

        return connected_components

    def _find_eulerian_path(self, points: Set[Point]) -> List[Point]:
        """Finds an Eulerian path through a set of points in the graph.

        An Eulerian path visits every edge exactly once. This method attempts
        to find such a path considering only edges between points in the input
        set.

        This function implements the Hierholzer's algorithm.

        See: https://en.wikipedia.org/wiki/Eulerian_path#Hierholzer's_algorithm

        Args:
            points: A set of points to find a path through.

        Returns:
            An ordered list of points representing the Eulerian path.

        Raises:
            NoEulerianPathError: If no Eulerian path exists through the given
                points.
        """
        outside = set(self._points.keys()) - points

        def get_neighbours(point: Point) -> Set[Point]:
            return self._points[point] - outside

        def count_edges(point: Point) -> int:
            return len(get_neighbours(point))

        def find_start_point() -> Point:
            odd_degree_points = [p for p in points if count_edges(p) % 2 == 1]
            if len(odd_degree_points) == 0:
                return next(iter(points))
            elif len(odd_degree_points) == 2:
                return odd_degree_points[0]
            else:
                print([(p, count_edges(p)) for p in points])
                raise NonEulerianPathError("No Eulerian path exists")

        remaining_edges = {p: get_neighbours(p).copy() for p in points}

        def find_path(start_point):
            """
            Moves from `start_point` until there is an edge that give us the
            possibility to move to another point. If there is no such edge,
            it returns the generated path.
            """
            current_path = [start_point]
            while len(remaining_edges[current_path[-1]]) > 0:
                next_point = remaining_edges[current_path[-1]].pop()
                remaining_edges[next_point].remove(current_path[-1])
                current_path.append(next_point)
            return current_path

        path = find_path(find_start_point())

        while any(remaining_edges.values()):
            for i, p in enumerate(path):
                if len(remaining_edges[p]) > 0:
                    adding_path = find_path(start_point=path[i])

                    assert adding_path[0] == path[i]
                    assert len(adding_path) > 1
                    assert adding_path[-1] == path[i]

                    path = path[:i] + adding_path + path[i + 1 :]
                    break

        return path

    def as_eulerian_paths(self) -> List[List[Point]]:
        """Convert the graph into a list of Eulerian paths.

        This method attempts to find Eulerian paths for each connected
        component in the graph.
        If the points of the graph are the vertices of the boundary of a region
        of a mask, such paths always exist and represent the boundary curves
        of the region.

        Returns:
            A list where each inner list represents an Eulerian path for a
                connected component. Each path is a sequence of points
                (x, y) that forms the path.

        Raises:
            NonEulerianPathError: If any component in the graph cannot be
                represented as an Eulerian path. This happens when the number of
                vertices with odd degree is not 0 or 2 in any component.
        """
        connected_components = self._get_connected_components()
        paths = []

        for component in connected_components:
            current_path = self._find_eulerian_path(component)
            paths.append(current_path)

        return paths


class ComponentMask2D(ComponentMask):
    """
    This is a specialized version of the ComponentMask class that adds some
    methods that are suitable only for the 2D case.

    The most relevant one is the `get_component_boundary_curves`, that returns
    the boundary of a component as a list of curves.
    """

    def __init__(self, mask: np.ndarray):
        if mask.ndim != 2:
            raise ValueError(
                f"mask must be a 2D array; received shape: {mask.shape}"
            )
        super().__init__(mask)

    def get_component_boundary_curves(
        self, component_id: int, closed_boundary=False
    ) -> List[List[Point]]:
        """Extracts the boundary curves of a component as lists of points.

        For a given component, this method returns its boundary as a collection
        of curves. Each curve is described by a sequence of vertex points that
        form a boundary polygon. Points are represented as integer coordinates
        (x, y). For a cell at position (i, j), its vertices are located at:
            - (i, j)     - top-left corner
            - (i + 1, j) - bottom-left corner
            - (i, j + 1) - top-right corner
            - (i + 1, j + 1) - bottom-right corner

        The number of returned curves depends on the component's topology:
        - A simple component (no holes, not touching domain boundary) returns
          one curve representing its external boundary
        - A component touching the domain boundary may return multiple curves,
          as boundary segments at domain edges are excluded by default
          (unless `closed_boundary=True`)
        - Components with holes return additional curves, one for each hole's
          boundary

        Args:
            component_id: Identifier of the target component for boundary
                computation
            closed_boundary: When True, includes points at domain edges. When
                False, only includes points between interior and exterior cells.
                Defaults to False.

        Returns:
            A list of boundary curves, where each curve is a list of (x, y) points
            forming a closed path. All coordinates are integers representing vertex
            positions.
        """
        component = self.get_component(component_id)
        boundary = self.get_component_boundary(
            component_id, closed_boundary=closed_boundary
        )
        boundary_cells = np.nonzero(boundary)

        graph = _PointGraph()

        # We fill the graph by adding an edge for each face;
        # For each cell of the component
        for cell_coords in zip(boundary_cells[0], boundary_cells[1]):
            # For each face of the cell; here we describe the face by identify
            # the neighbor cell in its direction. So we now that starting from
            # the indices `c = [i, j]`, the neighbor coordinates can be found
            # by imposing `c[axis] += side`
            for axis, side in cart_prod((0, 1), (-1, +1)):
                # Check if the current face is on the boundary of the domain
                on_boundary = False
                if side == -1 and cell_coords[axis] == 0:
                    on_boundary = True
                if side == +1 and cell_coords[axis] == boundary.shape[axis] - 1:
                    on_boundary = True

                # If we are not on the boundary of the domain, but we do not
                # have to add the faces on the boundary, we can skip this one
                if not closed_boundary and on_boundary:
                    continue

                # If this face is between two cells, we need to check if the
                # neighbor cell is outside the current component. If this is
                # the case, then the face is part of the boundary (and we need
                # to add it to our graph). Otherwise, we can just skip this
                # face.
                if not on_boundary:
                    neighbour_cell = list(cell_coords)
                    neighbour_cell[axis] += side
                    neighbour_cell = tuple(neighbour_cell)
                    if component[neighbour_cell]:
                        continue

                # This is a face that must be registered. We identify
                # the vertices of the face. `moving_axis` is the axis
                # orthogonal to the original `axis`, and it is the one along
                # which we move to find the other vertex.
                moving_axis = 1 - axis

                # In our first guess, the vertces of the face have the same
                # indices of the original cell. Moreover, we know that the
                # second point will have a +1 along the `moving_axis` axis.
                p1 = list(cell_coords)
                p2 = list(cell_coords)
                p2[moving_axis] += 1

                # If we are on the `right` or on the `top` of the cell, we must
                # also add +1 to the other axis.
                if side == 1:
                    p1[axis] += 1
                    p2[axis] += 1

                # Finally, we register the face in the graph. We explicitly
                # convert the coordinates to integers to avoid the numpy
                # types.
                graph.add_edge(
                    (int(p1[0]), int(p1[1])), (int(p2[0]), int(p2[1]))
                )

        return graph.as_eulerian_paths()
