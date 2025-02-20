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
