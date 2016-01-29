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
    def __init__(self, time_interval, archive_dir="/", postfix_dir="", glob_pattern="ave*gz"):
        """
        TimeSeries constructor.

        Args:
            - *time_interval*: a TimeInterval object that defines the interval
              of time you want to search in.
            - *archive_dir* (optional): path to the directory that contains the
              archived files (default: "/").
            - *postfix_dir* (optional): the last part of the path before the
              files that will be appended to the date (default: "").
              E.G.: "POSTPROC/AVE_FREQ_1".
            - *glob_pattern* (optional): glob-style pattern for the files
              (default: "ave*gz").
        """
        #Input validation
        if isinstance(time_interval, TimeInterval):
            self._time_interval = time_interval
        else:
            raise ValueError("time_interval must be an instance of TimeInterval")

        self._postfix_dir = str(postfix_dir)
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
        #Input validation
        if not isinstance(rundays, (list, tuple)):
            raise ValueError("rundays must be a list or a tuple")
        else:
            for e in rundays:
                if not isinstance(e, int) or (e < 1) or (e > 7):
                    raise ValueError("all the elements of rundays must be integers between 1 and 7.")
        output = list()
        t = self._time_interval.start_time
        while t <= self._time_interval.end_time:
            if t.isoweekday() in rundays:
                run_path = path.join(self._archive_dir, t.strftime('%Y%m%d'), self._postfix_dir)
                if path.exists(run_path) and path.isdir(run_path):
                    output.append((t, run_path))
            t += timedelta(1)
        return output

    def get_analysis_days(self, rundays=[2]):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2]).

        Returns: a list of tuples (datetime, filename) of assimilation/hindcast
        computations.
        """
        #Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays)
        output = list()
        #For each directory
        for directory in search_paths:
            #Get the files list
            file_list = glob(path.join(directory[1], self._glob_pattern))
            #Take the first seven days
            t = directory[0] - timedelta(7)
            #For each day
            for _ in range(7):
                #For each filename
                for filename in file_list:
                    #Get the basename
                    bn = path.basename(filename)
                    #If there's a filename with that date and within the time interval
                    if (bn.find(t.strftime("%Y%m%d")) != -1) and self._time_interval.contains(t):
                        #Build the tuple and append the tuple to the output
                        output.append((t,filename))
                #Get the next day
                t += timedelta(1)
        #Sort the output by date and return it to the caller
        output = sorted(output, key=lambda x: x[0])
        return output

    def get_forecast_days(self, rundays=[2,5]):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2,5]).

        Returns: a list of tuples (datetime, filename) of forecast computations.
        """
        #Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays)
        #Create the working dictionary
        wdict = dict()
        #For each directory
        for directory in search_paths:
            #Get the files list
            file_list = glob(path.join(directory[1], self._glob_pattern))
            #Start from the run day
            t = directory[0]
            #stop after 10 days
            stop_t = t + timedelta(10)
            #For each day
            while t < stop_t:
                #Build the datestring
                datestring = t.strftime("%Y%m%d")
                #For each file
                for filename in file_list:
                    #Get the base name
                    bn = path.basename(filename)
                    #If bn contains datestring
                    if bn.find(datestring) != -1:
                        #Add (t, filename) tuple to wdict
                        wdict[datestring] = (t, filename)
                #Increment t
                t += timedelta(1)
        #Sort wdict keys
        key_list = sorted(wdict.keys())
        #Create the output list
        output = list()
        #For each element of key_list
        for key in key_list:
            #Append wdict element to output
            output.append(wdict[key])
        return output
