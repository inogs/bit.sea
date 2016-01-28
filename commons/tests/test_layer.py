# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

from nose.tools import *
from ..layer import *

def test_initial():
    assert True

def test_layer_constructor():
    l = Layer(0, 10)
    assert not (l is None)
    assert (l.top == 0)
    assert (l.bottom == 10)

@raises(ValueError)
def test_layer_constructor_wrong_order():
    l = Layer(10, 0)

@raises(ValueError)
def test_layer_constructor_zero_height():
    l = Layer(10, 10)

def test_layer_repr():
    l = Layer(0, 10)
    assert (repr(l) == "Layer 0-10 m")

def test_layer_str():
    l = Layer(0, 10)
    assert (str(l) == "layer0-10")

def test_layer_string():
    l = Layer(0, 10)
    assert (l.string() == "0-10m")

def test_layer_longname():
    l = Layer(0, 10)
    assert (l.longname() == "0000-0010m")
