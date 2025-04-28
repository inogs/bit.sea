# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import subprocess
from datetime import datetime
from datetime import timedelta
from pathlib import Path
from string import Template
from warnings import warn

from bitsea.commons.time_interval import TimeInterval


class TimeSeries:
    """
    Time series handling class.
    """

    def __init__(
        self,
        time_interval,
        archive_dir="/",
        postfix_dir="",
        glob_pattern="ave*gz",
    ):
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
        # Input validation
        if isinstance(time_interval, TimeInterval):
            self._time_interval = time_interval
        else:
            raise ValueError(
                "time_interval must be an instance of TimeInterval"
            )

        self._postfix_dir = Path(postfix_dir)
        self._glob_pattern = str(glob_pattern)

        if self._postfix_dir.is_absolute():
            raise ValueError("postfix_dir must be a relative path")

        archive_dir = Path(archive_dir)
        if archive_dir.exists():
            if not archive_dir.is_dir():
                raise ValueError("%s is not a directory" % archive_dir)
            self._archive_dir = archive_dir.absolute()
        else:
            warn("%s doesn't exists" % archive_dir)
            self._archive_dir = archive_dir

    def get_runs(self, rundays=(2, 5), fudge=7):
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
        # Input validation
        if not isinstance(rundays, (list, tuple)):
            raise ValueError("rundays must be a list or a tuple")

        for e in rundays:
            if e < 1 or e > 7:
                raise ValueError(
                    "all the elements of rundays must be integers between 1"
                    "and 7."
                )

        output = list()

        t = self._time_interval.start_time
        stop_day = self._time_interval.end_time + timedelta(days=fudge)

        candidate_dates = []
        while t <= stop_day:
            candidate_dates.append(t)
            t += timedelta(days=1)

        for t in candidate_dates:
            if t.isoweekday() not in rundays:
                continue

            file_name = Path(t.strftime("%Y%m%d")) / self._postfix_dir
            run_path = self._archive_dir / file_name
            if run_path.is_dir():
                output.append((t, run_path))

        return output

    def get_analysis_days(self, rundays=(2,)):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2]).

        Returns: a list of tuples (datetime, filename) of assimilation/hindcast
        computations.
        """
        # Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays)
        output = list()
        # For each directory
        for time_obj, directory in search_paths:
            # Get the files list
            file_list = sorted(directory.glob(self._glob_pattern))

            # Take the first seven days
            farther_date = time_obj - timedelta(days=8)
            candidate_dates = tuple(
                farther_date + timedelta(days=i) for i in range(7)
            )

            for t in candidate_dates:
                # For each filename
                for file_path in file_list:
                    # Get the basename
                    bn = file_path.name
                    # If there's a filename with that date and within the time interval
                    if bn.find(
                        t.strftime("%Y%m%d")
                    ) != -1 and self._time_interval.contains(t):
                        # Build the tuple and append the tuple to the output
                        output.append((t, file_path))

        # Sort the output by date and return it to the caller
        output = sorted(output, key=lambda x: x[0])
        return output

    def get_simulation_days(self, rundays=(2,)):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2]).

        Returns: a list of tuples (datetime, filename) of assimilation/hindcast
        computations.
        """
        # Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays)
        output = list()
        assert self._postfix_dir != Path("OPAOPER_F")
        is_phys = self._postfix_dir == Path("OPAOPER_A")
        for t, directory in search_paths:
            # Get the files list
            file_list = directory.glob(self._glob_pattern)
            # For each day

            for file_path in file_list:
                # Get the basename
                bn = file_path.name
                # If there's a filename with that date and within the time interval
                if is_phys:
                    if (
                        bn.find(t.strftime("%Y%m%d_s_")) != -1
                    ) and self._time_interval.contains(t):
                        # Build the tuple and append the tuple to the output
                        output.append((t, file_path))
                else:
                    if (
                        bn.find(t.strftime("ave.%Y%m%d")) != -1
                    ) and self._time_interval.contains(t):
                        # Build the tuple and append the tuple to the output
                        output.append((t, file_path))

        # Sort the output by date and return it to the caller
        output = sorted(output, key=lambda x: x[0])
        return output

    def get_forecast_days(self, rundays=(2, 5)):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2,5]).

        Returns: a list of tuples (datetime, filename) of forecast computations.
        """
        is_phys = self._postfix_dir in (Path("OPAOPER_A"), Path("OPAOPER_F"))
        # Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays, fudge=0)

        # Create the working dictionary
        wdict = dict()

        # For each directory
        for start_time, directory in search_paths:
            local_dict = dict()
            # Get the files list
            file_list = sorted(directory.glob(self._glob_pattern))

            dates = tuple(start_time + timedelta(days=i) for i in range(10))

            for t in dates:
                for file_path in file_list:
                    # Get the base name
                    bn = file_path.name
                    date_string = t.strftime("%Y%m%d")
                    if is_phys:
                        if (
                            bn.find(t.strftime("%Y%m%d_f_")) != -1
                        ) and self._time_interval.contains(t):
                            # Build the tuple and append the tuple to the outputs
                            if date_string in local_dict:
                                local_dict[date_string].append(file_path)
                            else:
                                local_dict[date_string] = [file_path]
                    else:
                        if (
                            bn.find(t.strftime("%Y%m%d")) != -1
                        ) and self._time_interval.contains(t):
                            local_dict[date_string] = [file_path]

            wdict[directory] = local_dict

        unique_dict = dict()
        for t, directory in search_paths:
            local_dict = wdict[directory]
            for key in local_dict.keys():
                unique_dict[key] = local_dict[key]  # this can overwrite

        output = list()
        key_list = sorted(unique_dict.keys())
        for t in key_list:
            filelist = unique_dict[t]
            for file_path in filelist:
                output.append((t, file_path))

        return output

    def extract_analysis(
        self,
        outputdir,
        rundays=(2,),
        command="gzip -cd $INFILE > $OUTFILE",
        remove_ext=True,
    ):
        """
        Extracts analysis files to outputdir.

        Args:
            - *outputdir*: path to the output directory. If it does not exist
              this method will attempt to create it. If it exists, but it is not
              a directory a ValueError will be raised.
            - *rundays* (optional): see get_runs (default: [2]).
            - *command* (optional): the command string template. This template
              will be used to build the command string that will be passed to
              the shell. It must contain two literal substrings ($INFILE and
              $OUTFILE) that will be substituted with the input and output file
              paths (default: "gzip -cd $INFILE > $OUTFILE").
              E.G.:
                    cp $INFILE $OUTFILE
                    cat $INFILE > $OUTFILE
                    gzip -cd $INFILE > $OUTFILE
            - *remove_ext* (optional): boolean value. If set to True the output
              file basename will be equal to the input file base name without
              the extension, otherwise the base names will be the same
              (default: True).
              E.G. if remove_ext is True 'file01.gz' becomes 'file01'.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        file_list = self.get_analysis_days(rundays)
        outputdir = Path(outputdir).absolute()
        return self._extract(file_list, outputdir, str(command), remove_ext)

    def extract_forecast(
        self,
        outputdir,
        rundays=(2, 5),
        command="gzip -cd $INFILE > $OUTFILE",
        remove_ext=True,
    ):
        """
        Extracts forecast files to outputdir.

        Args: see extract_analysis.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        file_list = self.get_forecast_days(rundays)
        outputdir = Path(outputdir).absolute()
        return self._extract(file_list, outputdir, str(command), remove_ext)

    def extract_simulation(
        self,
        outputdir,
        rundays=(2,),
        command="gzip -cd $INFILE > $OUTFILE",
        remove_ext=True,
    ):
        """
        Extracts forecast files to outputdir.

        Args: see extract_analysis.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        file_list = self.get_simulation_days(rundays)
        outputdir = Path(outputdir).absolute()
        return self._extract(file_list, outputdir, str(command), remove_ext)

    def extract_from_list(
        self,
        file_list,
        outputdir,
        command="gzip -cd $INFILE > $OUTFILE",
        remove_ext=True,
    ):
        """
        Extract the files from a list into outputdir.

        Args:
            - *file_list*: a list of tuples (datetime, filename).
            - *outputdir*: path to the output directory. If it does not exist
              this method will attempt to create it. If it exists, but it is not
              a directory, a ValueError will be raised.
            - *command* (optional): the command string template. See extract_analysis.
            - *remove_ext* (optional): boolean value. See extract_analysis.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        return self._extract(file_list, outputdir, str(command), remove_ext)

    # Private methods
    def _extract(self, file_list, outputdir, command, remove_ext):
        outputdir = Path(outputdir)

        # If outputdir exists make sure it's a directory
        if outputdir.exists():
            if not outputdir.is_dir():
                raise ValueError("%s is not a directory" % outputdir)
        else:
            # Try to create the output directory
            outputdir.mkdir()

        # Build command template
        if (command.find("$INFILE") != -1) and (command.find("$OUTFILE") != -1):
            template = Template(command)
        else:
            raise ValueError("Invalid command string: " + command)

        output = list()
        # For each input file
        for _, ifn in file_list:
            # Get the path for the output file
            if remove_ext:
                file_name = ifn.stem
            else:
                file_name = ifn.name
            ofn = outputdir / file_name
            # Build the actual command string
            comstring = template.substitute(INFILE=ifn, OUTFILE=ofn)
            # Lauch the command
            retcode = subprocess.call(comstring, shell=True)
            # If the subprocess ended normally
            if retcode == 0:
                output.append(ofn)
            else:
                warn("Extract command exit code: %d" % retcode)
        return output

    # Static Methods
    @staticmethod
    def get_sublist(L, weekdays):
        """
        Filters a list by week day.

        Args:
            - *L*: a list of tuples (datetime, path).
            - *weekdays*: a list of numbers mapping the week days to take.
               1 Monday
               2 Tuesday
               3 Wednesday
               4 Thursday
               5 Friday
               6 Saturday
               7 Sunday

        Returns: a list of tuples (datetime, path) filtered according to
        weekdays. E.G.: if weekdays is [1,4] the resulting list will contain
        only Mondays and Thursdays.
        """
        if not isinstance(weekdays, (list, tuple)):
            raise ValueError("weekdays must be a list or a tuple")
        for el in weekdays:
            if not isinstance(el, int) or (el < 1) or (el > 7):
                raise ValueError("Invalid weekday: %s" % str(el))
        output = list()
        for el in L:
            if isinstance(el, tuple) and isinstance(el[0], datetime):
                if el[0].isoweekday() in weekdays:
                    output.append(el)
            else:
                warn("Skipping invalid element: %s" % str(el))
        return output

    @staticmethod
    def get_daily_sat(L, sat_archive, dateformat="%Y%m%d", glob_pattern="*.nc"):
        """
        Matches by date and file name the elements of L with the files into
        sat_archive directory.

        Args:
            - *L*: a list of tuples (datetime, path).
            - *sat_archive*: the path to the archive directory.
            - *dateformat* (optional): a string defining the date format (see
              datetime.strftime) to seek in the file name.

        Returns: a list of tuples (datetime, path, sat_file_path).
        """
        if not isinstance(L, (list, tuple)):
            raise ValueError("L must be a list or a tuple")
        sat_archive = Path(sat_archive)
        if sat_archive.exists():
            if not sat_archive.is_dir():
                raise ValueError("%s is not a directory" % sat_archive)
        else:
            raise ValueError("%s doesn't exist" % sat_archive)
        sat_files = sorted(sat_archive.glob(glob_pattern))
        output = list()
        for dt, p in L:
            ds = dt.strftime("%Y%m%d")
            for f in sat_files:
                bn = f.name
                if bn.find(ds) != -1:
                    output.append((dt, p, f))
        return output
