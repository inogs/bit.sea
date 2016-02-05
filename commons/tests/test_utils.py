# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import numpy as np

from nose.tools import *
from ..utils import *

def test_initial():
    assert True

def test_is_number():
    #Integer
    assert is_number(0)
    #Float
    assert is_number(0.0)
    #Complex
    assert is_number(0.0+1.0j)
    #Numpy float32
    assert is_number(np.float32(0))
    #None
    assert not is_number(None)
    #String
    assert not is_number("0")

def test_get_date_string():
    l, s = get_date_string("20160101")
    assert (l == "2016-01-01") and (s == "20160101")

def test_get_date_string_empty_string():
    l, s = get_date_string("")
    assert (l == "") and (l == s)

@raises(TypeError)
def test_get_date_string_invalid_input():
    #List of one string
        _ , _ = get_date_string(['20160101'])
    #Tuple with one string
        _ , _ = get_date_string(('20160101',))
    #Number
        _ , _ = get_date_string(20160101)
    #None
        _ , _ = get_date_string(None)
    #Empty list
        _ , _ = get_date_string([])
