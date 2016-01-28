# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from nose.tools import *
from ..timeseries import *
from ..time_interval import TimeInterval

def test_initial():
    assert True

def test_timeseries_constructor():
    ti = TimeInterval()
    assert not (ti is None)
    ts = TimeSeries(ti)
    assert not (ts is None)

@raises(ValueError)
def test_timeseries_constructor_wrong_ti():
    ts = TimeSeries(None)

@raises(ValueError)
def test_timeseries_constructor_wrong_archive():
    ti = TimeInterval()
    ts = TimeSeries(ti, archive_dir="/etc/hosts")

def test_timeseries_constructor_unexistant_archive():
    ti = TimeInterval()
    ts = TimeSeries(ti, archive_dir="doesnotexist")
    assert not (ts is None)

def test_timeseries_get_runs():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    o  = ts.get_runs()
    assert isinstance(o, list)

