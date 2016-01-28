# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import os
import shutil

from datetime import datetime, timedelta
from nose.tools import *
from ..timeseries import *
from ..time_interval import TimeInterval

def build_archive_files(directory):
    #Get the base date
    t = datetime.strptime(os.path.basename(directory), '%Y%m%d')
    t -= timedelta(7)
    #Stop date
    st = t + timedelta(17)
    while t < st:
        filename = os.path.join(directory, "%s%s%s" % ("ave", t.strftime('%Y%m%d'), ".gz"))
        fd = open(filename, 'w')
        fd.close()
        t += timedelta(1)

def build_tmp_archive():
    if not os.path.exists('/tmp/archive'):
        os.mkdir('/tmp/archive')
        os.mkdir('/tmp/archive/20160101')
        build_archive_files('/tmp/archive/20160101')
        os.mkdir('/tmp/archive/20160105')
        build_archive_files('/tmp/archive/20160105')
        os.mkdir('/tmp/archive/20160108')
        build_archive_files('/tmp/archive/20160108')

def clean_tmp_archive():
    if os.path.exists('/tmp/archive'):
        shutil.rmtree('/tmp/archive')

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

@raises(ValueError)
def test_timeseries_get_runs_rundays_not_a_list():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    o  = ts.get_runs(rundays=None)

@raises(ValueError)
def test_timeseries_get_runs_rundays_number():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    o  = ts.get_runs(rundays=2)

@raises(ValueError)
def test_timeseries_get_runs_rundays_overflow():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    o  = ts.get_runs(rundays=[10,20])

def test_timeseries_get_analysis_days():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.get_analysis_days()
    assert isinstance(o, list)
    assert len(o) > 0
    clean_tmp_archive()

def test_timeseries_get_forecast_days():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.get_forecast_days()
    assert isinstance(o, list)
    assert len(o) > 0
    clean_tmp_archive()
