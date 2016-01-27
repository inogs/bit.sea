# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

from os import path
from glob import glob
from datetime import datetime, timedelta
from warnings import warn

from commons.time_interval import TimeInterval

class TimeSeries(object):
    """
    Time series handling class.
    """
    def __init__(self, time_interval, archive_dir="/", prefix_dir="", glob_pattern="ave*gz"):
        #Input validation
        if isinstance(time_interval, TimeInterval):
            self._time_interval = time_interval
        else:
            raise ValueError("time_interval must be an instance of TimeInterval")

        self._prefix_dir = str(prefix_dir)
        self._glob_pattern = str(glob_pattern)

        archive_dir = str(archive_dir)
        if path.exists(archive_dir):
            if path.isdir(archive_dir):
                self._archive_dir = path.abspath(archive_dir)
            else:
                raise ValueError("%s is not a directory" % archive_dir)
        else:
            warn("%s doesn't exists" % archive_dir)
            self._archive_dir = archive_dir

    def get_runs(self, rundays=[2,5]):
        """
        Args:
            - *rundays* (optional): a list of integers mapping the weekdays in
              which the chain runs (default: [2,5]).
              1 is Monday
              2 is Tuesday
              3 is Wednesday
              4 is Thursday
              5 is Friday
              6 is Saturday
              7 is Sunday

        Returns: a list of tuples (datetime, path) within the time_interval for
        which there's a run.
        """
        output = list()
        t = self._time_interval.start_time
        while t <= self._time_interval.end_time:
            if t.isoweekday() in rundays:
                run_path = path.join(self._archive_dir, t.strftime('%Y%m%d'), self._prefix_dir)
                if path.exists(run_path) and path.isdir(run_path):
                    output.append((t, run_path))
            t += timedelta(1)
        return output
