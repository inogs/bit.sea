from nose.tools import *
from ..mean import Mean as m
from ..mean import GaussianMean as gm
from numpy import array
from numpy import linspace
from numpy import nan
from numpy import isnan

def test_initial():
    assert True

#mean base class tests
@raises(NotImplementedError)
def test_create_mean_object():
    a = m(0)

@raises(NotImplementedError)
def test_mean_object_compute():
    a = m(0)
    a.compute([1,])

#Gaussian mean tests
def test_create_gmean_object():
    a = gm(0)

@raises(ValueError)
def test_create_gmean_object_with_negative_sigma():
    a = gm(0, -1.0)

@raises(ValueError)
def test_create_gmean_object_with_sigma_string():
    a = gm(0, 'hello')

@raises(ValueError)
def test_create_gmean_object_negative():
    a = gm(-1)

@raises(ValueError)
def test_create_gmean_object_with_string():
    a = gm('hello')

def test_create_gmean_object_with_float():
    a = gm(1.0)

@raises(TypeError)
def test_gm_compute_no_params():
    a = gm(0)
    a.compute()

@raises(TypeError)
def test_gm_compute_single_digits():
    a = gm(0)
    a.compute(1)

@raises(TypeError)
def test_gm_compute_strings():
    a = gm(0)
    a.compute('1')

def test_gm_compute_lists_of_single_point():
    a = gm(0)
    r = a.compute([1,])
    assert len(r) == 1
    assert r[0] == 1

def test_gm_compute_tuples_of_single_point():
    a = gm(0)
    r = a.compute((1,))
    assert len(r) == 1
    assert r[0] == 1

def test_gm_compute_nparrays_of_single_point():
    obj = gm(0)
    a = array([1,])
    result = obj.compute(a)
    assert len(result) == 1
    assert result[0] == 1

def test_gm_compute_nparrays_one_point_interval():
    obj = gm(0)
    a = array([1, 2, 3])
    result = obj.compute(a)
    assert len(result) == len(a)
    for i in range(len(a)):
        assert result[i] == a[i]

def test_gm_compute_nparrays_three_points_interval():
    obj = gm(3)
    a = array([1, 1, 1])
    result = obj.compute(a)
    assert len(result) == len(a)
    for i in range(len(a)):
        assert abs(result[i] - a[i]) < 1e-16

#Test numerical stability
def test_gm_compute_1001nparray_symmetrical():
    #Generate symmetric array of 1001 elements from -500 to 500 included
    #a[500] is always 0
    a = linspace(-500, 500, 1001)
    assert a[500] == 0.0
    #Test 3,5,7,9,11,13 and 15 elements-long intervals
    for i in range(3,16,2):
        obj = gm(i)
        result = obj.compute(a)
        assert len(result) == len(a)
        assert abs(result[500]) < 1e-16
