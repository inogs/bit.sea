# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

from nose.tools import *

from ..segment import *

def test_initial():
    assert True

def test_segment_constructor():
    s = Segment((0.0,0.0), (0.0,90.0))
    assert not (s is None)
    assert (s.lon_min == 0.0)
    assert (s.lon_max == 0.0)
    assert (s.lat_min == 0.0)
    assert (s.lat_max == 90.0)
    assert (s.lon_lat_min == (0.0,0.0))
    assert (s.lon_lat_max == (0.0,90.0))

@raises(ValueError)
def test_segment_constructor_none():
    s = Segment(None, None)

@raises(ValueError)
def test_segment_constructor_one_element_per_tuple():
    s = Segment((0.0,), (0.0,))

@raises(ValueError)
def test_segment_constructor_nan():
    s = Segment((0.0,None), (0.0,None))

def test_segment_constructor_diagonal():
    s = Segment((0.0,0.0), (45.0,45.0))
    assert not (s is None)

@raises(ValueError)
def test_segment_constructor_diagonal_zero_points():
    s = Segment((0.0,0.0), (45.0,45.0), points=0)

def test_segment_name():
    s = Segment((0.0,0.0),(0.0,90.0), name="A")
    assert (s.name == "A")

def test_segment_str():
    s = Segment((0.0,0.0), (0.0,90.0))
    assert (str(s) == "Segment from 0,0 to 0,90")
    s = Segment((0.0,0.0),(0.0,90.0), name="A")
    assert (str(s) == "Segment A from 0,0 to 0,90")

def test_segment_repr():
    s = Segment((0.0,0.0), (0.0,90.0))
    assert (repr(s) == "Segment((0.000000,0.000000), (0.000000,90.000000))")
    s = Segment((0.0,0.0),(0.0,90.0), name="A")
    assert (repr(s) == "Segment A ((0.000000,0.000000), (0.000000,90.000000))")

