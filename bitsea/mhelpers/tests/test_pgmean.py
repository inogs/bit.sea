# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from nose.tools import *
from ..pgmean import PGaussianMean as gm
from ..pgmean import PLGaussianMean as pgm
from numpy import array
from numpy import linspace
from numpy import nan
from numpy import isnan

epsilon = 2.22e-16

def test_initial():
    assert True

def test_create_pgmean_object():
    obj = gm(0)

def test_pgm_compute_one_element_list():
    obj = gm(3)
    a = [1]
    b = [1]
    result = obj.compute(a, b)
    assert len(result) == len(a)
    assert result[0] == a[0]

def test_pgm_compute_nparrays_three_points_interval():
    obj = gm(3)
    a = array([1, 1, 1])
    b = array([1, 2, 3])
    result = obj.compute(a, b)
    assert len(result) == len(a)
    for i in range(len(a)):
        assert abs(result[i] - a[i]) < epsilon

@raises(AssertionError)
def test_pgm_compute_nparray_None():
    obj = gm(3)
    a = array([1, 1, 1])
    result = obj.compute(a, None)

#Test numerical stability
def test_pgm_compute_1001nparray_symmetrical():
    #Generate symmetric array of 1001 elements from -500 to 500 included
    #a[500] is always 0
    a = linspace(-500, 500, 1001)
    assert a[500] == 0.0
    #Generate pressures array
    b = linspace(0,1001,1001)
    #Test 3,5,7,9,11,13 and 15 elements-long intervals
    for i in range(3,16,2):
        obj = gm(i)
        result = obj.compute(a, b)
        assert len(result) == len(a)
        assert abs(result[500]) < epsilon

#PLGaussianMean Tests
def test_create_plgmean_object():
    obj = pgm(0)

@raises(AssertionError)
def test_plgm_compute_null_pressures():
    obj = pgm(3)
    obj.compute([1,2,3], None)

def test_plgm_compute_1001nparray_symmetrical():
    epsilon = 1e-14
    #Generate symmetric array of 1001 elements from -500 to 500 included
    #a[500] is always 0
    a = linspace(-500, 500, 1001)
    assert a[500] == 0.0
    #Generate pressures array
    b = linspace(0,1001,1001)
    #Test 3,5,7,9,11,13 and 15 elements-long intervals
    for i in range(3,16,2):
        obj = pgm(i)
        result = obj.compute(a, b)
        assert len(result) == len(a)
        if abs(result[500]) > epsilon:
            raise ValueError("result (%s) is greater than epsilon (%e) by %e" % (abs(result[500]), epsilon, (abs(result[500]) - epsilon)))

