# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import os
import shutil
import subprocess

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
        p = subprocess.Popen(['gzip','-c'], stdin=subprocess.PIPE, stdout=fd)
        p.communicate(t.strftime('%Y-%m-%d'))
        p.wait()
        fd.close()
        t += timedelta(1)

def build_sat_archive(directory, timeinterval):
    if not os.path.exists(directory):
        os.mkdir(directory)
    t = timeinterval.start_time
    while t < timeinterval.end_time:
        filename = os.path.join(directory, "%s%s%s" % ("sat", t.strftime('%Y%m%d'), ".nc"))
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

def clean_tmp_archive(tmpdir='/tmp/archive'):
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)

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
    assert len(o) == 4
    clean_tmp_archive()

def test_timeseries_get_forecast_days():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.get_forecast_days()
    assert isinstance(o, list)
    assert len(o) == 17
    clean_tmp_archive()

def test_timeseries_extract_analysis():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.extract_analysis('/tmp/output')
    assert isinstance(o, list)
    assert len(o) == 4
    clean_tmp_archive()
    clean_tmp_archive('/tmp/output')

def test_timeseries_extract_forecast():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.extract_forecast('/tmp/output')
    assert isinstance(o, list)
    assert len(o) == 17
    clean_tmp_archive()
    clean_tmp_archive('/tmp/output')

@raises(ValueError)
def test_timeseries_extract_wrong_outputdir():
    ts = TimeSeries(TimeInterval())
    ts._extract([], '/etc/hosts', '', True)

@raises(ValueError)
def test_timeseries_extract_wrong_command():
    ts = TimeSeries(TimeInterval())
    ts._extract([], '/tmp', '', True)

def test_timeseries_extract_retcode_nonzero():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.extract_forecast("/tmp/archive", command="echo $INFILE $OUTFILE && exit 2")
    clean_tmp_archive()

def test_timeseries_extract_preserve_ext():
    build_tmp_archive()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir="/tmp/archive")
    o = ts.extract_forecast("/tmp/output", command="cp $INFILE $OUTFILE", remove_ext=False)
    assert len(o) == 17
    clean_tmp_archive()
    clean_tmp_archive('/tmp/output')

def test_timeseries_get_sublist():
    L = list()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    t = ti.start_time
    while t < ti.end_time:
        L.append((t,"",""))
        t += timedelta(1)
    L.append(None)
    L.append([t,"false","false"])
    o = TimeSeries.get_sublist(L, [1])
    assert isinstance(o, list)
    assert len(o) == 4

@raises(ValueError)
def test_timeseries_get_sublist_weekdays_not_list():
    TimeSeries.get_sublist([], 1)

@raises(ValueError)
def test_timeseries_get_sublist_weekdays_wrong_list():
    TimeSeries.get_sublist([(datetime.now, "")], [None])

def test_timeseries_get_daily_sat():
    L = list()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    t = ti.start_time
    while t < ti.end_time:
        L.append((t,"file_%s" % t.strftime("%Y%m%d")))
        t += timedelta(1)
    build_sat_archive("/tmp/sat_archive", ti)
    OL = TimeSeries.get_daily_sat(L, "/tmp/sat_archive")
    assert isinstance(OL, list)
    assert len(L) == len(OL)
    assert OL[0] == (L[0][0], L[0][1], "/tmp/sat_archive/sat%s.nc" % L[0][0].strftime("%Y%m%d"))
    clean_tmp_archive("/tmp/sat_archive")

@raises(ValueError)
def test_timeseries_get_daily_sat_L_not_list():
    TimeSeries.get_daily_sat(None, "")

@raises(ValueError)
def test_timeseries_get_daily_sat_sat_archive_not_string():
    TimeSeries.get_daily_sat([], None)

@raises(ValueError)
def test_timeseries_get_daily_sat_sat_archive_not_dir():
    TimeSeries.get_daily_sat([], "/etc/hosts")

def test_timeseries_get_daily_sat_empty_list():
    OL = TimeSeries.get_daily_sat([], "/")
    assert isinstance(OL, list)
    assert len(OL) == 0
