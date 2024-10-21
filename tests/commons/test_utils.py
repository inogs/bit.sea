# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import pytest

from bitsea.commons.utils import get_date_string


def test_get_date_string():
    l, s = get_date_string("20160101")
    assert (l == "2016-01-01") and (s == "20160101")

def test_get_date_string_empty_string():
    l, s = get_date_string("")
    assert (l == "") and (l == s)
