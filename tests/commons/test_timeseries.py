# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import gzip
import pytest
import shutil
from datetime import datetime
from datetime import date
from datetime import timedelta
from pathlib import Path

from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.timeseries import TimeSeries


def build_archive_files(dir_path: Path):
    dir_path.mkdir(exist_ok=True)
    #Get the base date
    dir_datetime = datetime.strptime(dir_path.name, '%Y%m%d')
    dir_date = date(dir_datetime.year, dir_datetime.month, dir_datetime.day)

    start_date = dir_date - timedelta(days=7)
    dates = tuple(
        start_date + timedelta(days=i) for i in range(17)
    )

    for t in dates:
        file_name = f"ave{t.strftime('%Y%m%d')}.gz"
        file_path = dir_path / file_name

        with gzip.open(file_path, 'wt') as fd:
            fd.write(t.strftime('%Y-%m-%d'))


def build_sat_archive(directory: Path, timeinterval):
    directory.mkdir(exist_ok=True)
    t = timeinterval.start_time
    while t < timeinterval.end_time:
        file_name = f"sat{t.strftime('%Y%m%d')}.nc"
        file_path = directory / file_name
        file_path.touch()
        t += timedelta(1)


@pytest.fixture(scope='module')
def test_archives(tmpdir_factory):
    archive_dir = Path(tmpdir_factory.mktemp("archive"))

    for dates in ("20160101", "20160105", "20160108"):
        build_archive_files(archive_dir / dates)

    yield archive_dir

    shutil.rmtree(archive_dir)


def test_timeseries_constructor():
    ti = TimeInterval()
    assert ti is not None
    ts = TimeSeries(ti)
    assert ts is not None


def test_timeseries_constructor_wrong_ti():
    with pytest.raises(ValueError):
        ts = TimeSeries(None)


def test_timeseries_constructor_wrong_archive():
    ti = TimeInterval()
    with pytest.raises(ValueError):
        ts = TimeSeries(ti, archive_dir="/etc/hosts")


def test_timeseries_constructor_nonexistent_archive():
    ti = TimeInterval()
    with pytest.warns(UserWarning):
        ts = TimeSeries(ti, archive_dir="doesnotexist")
    assert ts is not None


def test_timeseries_get_runs():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    o  = ts.get_runs()
    assert isinstance(o, list)


def test_timeseries_get_runs_rundays_not_a_list():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    with pytest.raises(ValueError):
        ts.get_runs(rundays=None)


def test_timeseries_get_runs_rundays_number():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    with pytest.raises(ValueError):
        ts.get_runs(rundays=2)


def test_timeseries_get_runs_rundays_overflow():
    ti = TimeInterval()
    ts = TimeSeries(ti)
    with pytest.raises(ValueError):
        ts.get_runs(rundays=[10,20])


def test_timeseries_get_analysis_days(test_archives):
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir=test_archives)
    o = ts.get_analysis_days()
    assert isinstance(o, list)
    assert len(o) == 3


def test_timeseries_get_forecast_days(test_archives):
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir=test_archives)
    o = ts.get_forecast_days()
    assert isinstance(o, list)
    assert len(o) == 17


def test_timeseries_extract_analysis(test_archives, tmp_path):
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir=test_archives)
    o = ts.extract_analysis(tmp_path)

    assert isinstance(o, list)
    assert len(o) == 3


def test_timeseries_extract_forecast(test_archives, tmp_path):
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir=test_archives)
    o = ts.extract_forecast(tmp_path)
    assert isinstance(o, list)
    assert len(o) == 17


def test_timeseries_extract_wrong_outputdir():
    ts = TimeSeries(TimeInterval())
    with pytest.raises(ValueError):
        ts._extract([], '/etc/hosts', '', True)


def test_timeseries_extract_wrong_command():
    ts = TimeSeries(TimeInterval())
    with pytest.raises(ValueError):
        ts._extract([], '/tmp', '', True)


def test_timeseries_extract_retcode_nonzero(test_archives):
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir=test_archives)
    with pytest.warns(UserWarning):
        ts.extract_forecast(
            test_archives,
            command="echo $INFILE $OUTFILE && exit 2"
        )

def test_timeseries_extract_preserve_ext(test_archives, tmp_path):
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    ts = TimeSeries(ti, archive_dir=test_archives)
    o = ts.extract_forecast(tmp_path, command="cp $INFILE $OUTFILE", remove_ext=False)
    assert len(o) == 17


def test_timeseries_get_sublist():
    L = list()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    t = ti.start_time
    while t < ti.end_time:
        L.append((t,"",""))
        t += timedelta(1)

    L.append(None)
    L.append([t, "false", "false"])

    with pytest.warns(UserWarning):
        o = TimeSeries.get_sublist(L, [1])
    assert isinstance(o, list)
    assert len(o) == 4


def test_timeseries_get_sublist_weekdays_not_list():
    with pytest.raises(ValueError):
        TimeSeries.get_sublist([], 1)


def test_timeseries_get_sublist_weekdays_wrong_list():
    with pytest.raises(ValueError):
        TimeSeries.get_sublist([(datetime.now, "")], [None])


def test_timeseries_get_daily_sat(tmp_path):
    time_list = list()
    ti = TimeInterval(starttime="20160101", endtime="20160131")
    t = ti.start_time
    while t < ti.end_time:
        time_list.append((t, f"file_{t.strftime('%Y%m%d')}"))
        t += timedelta(1)
    build_sat_archive(tmp_path, ti)
    OL = TimeSeries.get_daily_sat(time_list, tmp_path)
    assert isinstance(OL, list)
    assert len(time_list) == len(OL)

    expected_file = tmp_path / f'sat{time_list[0][0].strftime("%Y%m%d")}.nc'
    assert OL[0] == (time_list[0][0], time_list[0][1], expected_file)


def test_timeseries_get_daily_sat_L_not_list():
    with pytest.raises(ValueError):
        TimeSeries.get_daily_sat(None, "")


def test_timeseries_get_daily_sat_sat_archive_not_string():
    with pytest.raises(TypeError):
        TimeSeries.get_daily_sat([], None)


def test_timeseries_get_daily_sat_sat_archive_not_dir():
    with pytest.raises(ValueError):
        TimeSeries.get_daily_sat([], "/etc/hosts")


def test_timeseries_get_daily_sat_empty_list():
    daily_sat = TimeSeries.get_daily_sat([], "/")
    assert isinstance(daily_sat, list)
    assert len(daily_sat) == 0
