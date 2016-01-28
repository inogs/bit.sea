# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import numpy as np

from nose.tools import *
from ..helpers import *

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
