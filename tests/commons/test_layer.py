# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import pytest

from bitsea.commons.layer import Layer


def test_layer_constructor():
    layer = Layer(0, 10)
    assert layer.top == 0
    assert layer.bottom == 10


def test_layer_constructor_wrong_order():
    with pytest.raises(ValueError):
        Layer(10, 0)


def test_layer_constructor_zero_height():
    with pytest.raises(ValueError):
        Layer(10.1, 10)


def test_layer_repr():
    layer = Layer(0, 10)
    assert layer == eval(repr(layer))


def test_layer_str():
    layer = Layer(0, 10)
    assert str(layer) == "Layer 0-10 m"


def test_layer_string():
    layer = Layer(0, 10)
    assert layer.string() == "0-10m"


def test_layer_longname():
    layer = Layer(0, 10)
    assert layer.longname() == "0000-0010m"
