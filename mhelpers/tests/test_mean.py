from nose.tools import *
from .. import mean as m
from numpy import array
from numpy import nan
from numpy import isnan

def test_initial():
    assert True

#Test if you can create a mean object
def test_create_mean_object():
    a = m.mean(0)

@raises(ValueError)
def test_create_mean_object_with_negative_sigma():
    a = m.mean(0, -1.0)

@raises(ValueError)
def test_create_mean_object_with_sigma_string():
    a = m.mean(0, 'hello')

@raises(ValueError)
def test_create_mean_object_negative():
    a = m.mean(-1)

@raises(ValueError)
def test_create_mean_object_with_string():
    a = m.mean('hello')

def test_create_mean_object_with_float():
    a = m.mean(1.0)

@raises(TypeError)
def test_compute_no_params():
    a = m.mean(0)
    a.compute()

@raises(TypeError)
def test_compute_single_digits():
    a = m.mean(0)
    a.compute(1)

@raises(TypeError)
def test_compute_strings():
    a = m.mean(0)
    a.compute('1')

def test_compute_lists_of_single_point():
    a = m.mean(0)
    r = a.compute([1,])
    assert len(r) == 1
    assert r[0] == 1

def test_compute_tuples_of_single_point():
    a = m.mean(0)
    r = a.compute((1,))
    assert len(r) == 1
    assert r[0] == 1

def test_compute_nparrays_of_single_point():
    obj = m.mean(0)
    a = array([1,])
    result = obj.compute(a)
    assert len(result) == 1
    assert result[0] == 1

def test_compute_nparrays_one_point_interval():
    obj = m.mean(0)
    a = array([1, 2, 3])
    result = obj.compute(a)
    assert len(result) == len(a)
    for i in range(len(a)):
        assert result[i] == a[i]

def test_compute_nparrays_three_points_interval():
    obj = m.mean(3)
    a = array([1, 1, 1])
    result = obj.compute(a)
    assert len(result) == len(a)
    for i in range(len(a)):
        assert abs(result[i] - a[i]) < 1e-16

