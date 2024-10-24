# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import pytest

from bitsea.commons.layer import Layer


def test_layer_constructor():
    l = Layer(0, 10)
    assert (l.top == 0)
    assert (l.bottom == 10)


def test_layer_constructor_wrong_order():
    with pytest.raises(ValueError):
        Layer(10, 0)

def test_layer_constructor_zero_height():
    with pytest.raises(ValueError):
        Layer(10.1, 10)

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
