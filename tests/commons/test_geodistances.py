import numpy as np

from bitsea.commons.geodistances import extend_from_average
from bitsea.commons.geodistances import GeoPySidesCalculator
from bitsea.commons.geodistances import NemoGridSidesCalculator


def test_different_algorithms_give_coherent_results():
    grid_shape = (75, 130)

    xlevels = np.linspace(-5, 5, grid_shape[1])
    xlevels = xlevels * np.linspace(1, 2, grid_shape[0])[:, np.newaxis]

    y_levels = np.linspace(20, 30, grid_shape[0])
    ylevels = np.broadcast_to(y_levels[:, np.newaxis], grid_shape)

    c1 = NemoGridSidesCalculator()
    c2 = GeoPySidesCalculator(geodesic=False)
    c3 = GeoPySidesCalculator(geodesic=True)

    e1t_c1, e2t_c1 = c1(xlevels, ylevels)
    e1t_c2, e2t_c2 = c2(xlevels, ylevels)
    e1t_c3, e2t_c3 = c3(xlevels, ylevels)

    assert np.allclose(e1t_c1, e1t_c2, rtol=1e-2)
    assert np.allclose(e1t_c2, e1t_c3, rtol=1e-2)

    assert np.allclose(e2t_c1, e2t_c2, rtol=1e-1)
    assert np.allclose(e2t_c2, e2t_c3, rtol=1e-2)


def test_extend_from_average():
    t_array = np.linspace(0, 10, 7)
    t_array *= t_array

    v_array = extend_from_average(t_array)
    v_mean = (v_array[1:] + v_array[:-1]) / 2

    assert np.allclose(v_mean, t_array)