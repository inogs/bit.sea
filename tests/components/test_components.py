import numpy as np

from bitsea.components.component_mask import ComponentMask
from bitsea.components.find_component import find_component


def test_find_components_2d():
    data_map = np.zeros((10, 10), dtype=bool)
    data_map[1:4, 1:4] = True
    data_map[5:8, 1:9] = True
    data_map[1:9, 5:9] = True

    component1 = find_component(data_map, (2, 2))
    component2 = find_component(data_map, (7, 7))
    assert not np.any(np.logical_and(component1, component2))
    assert np.all(np.logical_or(component1, component2) == data_map)


def test_find_components_2d_with_border():
    data_map = np.zeros((7, 7), dtype=bool)
    data_map[0:2, 0:2] = True
    data_map[0:2, 5:] = True
    data_map[5:, 0:2] = True
    data_map[5:, 5:] = True
    data_map[1, 1:3] = True
    data_map[1:6, 3] = True
    data_map[5, 3:6] = True

    component = find_component(data_map, (0, 0))
    right_corner = find_component(data_map, (0, -1))
    left_corner = find_component(data_map, (-1, 0))
    corners = np.logical_or(right_corner, left_corner)

    # The components are separated
    assert not np.any(np.logical_and(component, corners))

    # Corners and component are all the original map
    assert np.all(np.logical_or(component, corners) == data_map)


def test_find_components_3d():
    data_map = np.zeros((3, 10, 10), dtype=bool)
    data_map[:, 1:4, 1:4] = True
    data_map[:, 5:8, 1:9] = True
    data_map[:, 1:9, 5:9] = True
    data_map[1, 3:6, 3:6] = True

    component1 = find_component(data_map, (0, 2, 2))
    component2 = find_component(data_map, (2, 7, 7))
    assert np.all(component1 == component2)
    assert np.all(component1 == data_map)


def test_component_mask_2d():
    data_map = np.zeros((7, 7), dtype=bool)
    data_map[0:2, 0:2] = True
    data_map[0:2, 5:] = True
    data_map[5:, 0:2] = True
    data_map[5:, 5:] = True
    data_map[1, 1:3] = True
    data_map[1:6, 3] = True
    data_map[5, 3:6] = True

    component_mask = ComponentMask(data_map)
    assert component_mask.n_components == 3
    assert np.all(np.equal(component_mask.components >= 0, data_map))


def test_component_mask_3d():
    data_map = np.zeros((3, 10, 10), dtype=bool)
    data_map[:, 1:4, 1:4] = True
    data_map[:, 5:8, 1:9] = True
    data_map[:, 1:9, 5:9] = True
    data_map[1, 3:6, 3:6] = True

    component_mask = ComponentMask(data_map)
    assert component_mask.n_components == 1
    assert np.all(component_mask.get_component(0) == data_map)


def test_component_mask_compute_boundary():
    data_map = np.zeros((10, 10, 10), dtype=bool)
    data_map[7, 1:9, 1:9] = True
    data_map[1:4, 1:4, 1:4] = True

    component_mask = ComponentMask(data_map)
    assert component_mask.n_components == 2

    k1 = component_mask.get_biggest_component()
    k2 = 1 if k1 == 0 else 0

    b2 = np.zeros_like(data_map)
    b2[1, 1:4, 1:4] = True
    b2[3, 1:4, 1:4] = True
    b2[1:4, 1, 1:4] = True
    b2[1:4, 3, 1:4] = True
    b2[1:4, 1:4, 1] = True
    b2[1:4, 1:4, 3] = True

    assert np.all(
        component_mask.get_component_boundary(k1)
        == component_mask.get_component(k1)
    )
    assert np.all(component_mask.get_component_boundary(k2) == b2)


def test_component_mask_compute_boundary_on_boundary():
    data_map = np.zeros((6, 11), dtype=bool)
    data_map[2:5, 2:9] = True
    data_map[2:5, 2:9] = True
    data_map[:3, 4:7] = True

    expected_boundary = np.zeros_like(data_map)
    expected_boundary[2:5, 2] = True
    expected_boundary[2:5, 8] = True
    expected_boundary[4, 2:9] = True
    expected_boundary[2, 2:5] = True
    expected_boundary[2, 6:9] = True
    expected_boundary[:2, 4] = True
    expected_boundary[:2, 6] = True

    expected_boundary_on_boundary = expected_boundary.copy()
    expected_boundary_on_boundary[0, 5] = True

    component_mask = ComponentMask(data_map)

    assert np.all(component_mask.get_component_boundary(0) == expected_boundary)
    assert np.all(
        component_mask.get_component_boundary(0, closed_boundary=True)
        == expected_boundary_on_boundary
    )


def test_component_mask_compute_external_boundary():
    data_map = np.zeros((11, 11), dtype=bool)
    data_map[2:9, 2:9] = True
    data_map[3:5, 3:5] = False
    data_map[6:8, 6:8] = False
    data_map[2, 2] = False
    data_map[8, 6] = False

    expected_boundary_with_holes = data_map.copy()
    expected_boundary_with_holes[3:5, 6:8] = False
    expected_boundary_with_holes[6:8, 3:5] = False

    expected_external_boundary = data_map.copy()
    expected_external_boundary[3, 3] = True
    expected_external_boundary[5:8, 3:5] = False
    expected_external_boundary[3:5, 5:8] = False

    component_mask = ComponentMask(data_map)

    assert np.all(
        component_mask.get_component_boundary(0) == expected_boundary_with_holes
    )
    assert np.all(
        component_mask.get_component_external_boundary(0)
        == expected_external_boundary
    )
