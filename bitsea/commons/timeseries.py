# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import subprocess as SP

from os import path, mkdir
from glob import glob
from datetime import datetime, timedelta
from warnings import warn
from string import Template

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

    def get_runs(self, rundays=[2,5], fudge=7):
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
        stop_day = self._time_interval.end_time + timedelta(days=fudge)
        while t <= stop_day:
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
        for time_obj, directory in search_paths:
            #Get the files list
            file_list = glob(path.join(directory, self._glob_pattern))
            #Take the first seven days
            t = time_obj - timedelta(8)
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

    def get_simulation_days(self,rundays=[2]):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2]).

        Returns: a list of tuples (datetime, filename) of assimilation/hindcast
        computations.
        """
        #Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays)
        output = list()
        assert self._postfix_dir !="OPAOPER_F/"
        is_phys = self._postfix_dir=='OPAOPER_A/'
        for directory in search_paths:
            #Get the files list
            file_list = glob(path.join(directory[1], self._glob_pattern))
            #Take the first seven days
            t = directory[0]
            #For each day
            
            for filename in file_list:
                #Get the basename
                bn = path.basename(filename)
                #If there's a filename with that date and within the time interval
                if is_phys:
                    if (bn.find(t.strftime("%Y%m%d_s_")) != -1) and self._time_interval.contains(t):
                        #Build the tuple and append the tuple to the output
                        output.append((t,filename))
                else:
                    if (bn.find(t.strftime("ave.%Y%m%d")) != -1) and self._time_interval.contains(t):
                        #Build the tuple and append the tuple to the output
                        output.append((t,filename))
                    

        #Sort the output by date and return it to the caller
        output = sorted(output, key=lambda x: x[0])
        return output
        #For each directory        
    def get_forecast_days(self, rundays=[2,5]):
        """
        Args:
            - *rundays* (optional): see get_runs (default: [2,5]).

        Returns: a list of tuples (datetime, filename) of forecast computations.
        """
        is_phys = (self._postfix_dir =="OPAOPER_A/") or (self._postfix_dir=='OPAOPER_F/')
        #Build the list of paths where we have to search for the files
        search_paths = self.get_runs(rundays, fudge=0)
        #Create the working dictionary
        wdict = dict()
        #For each directory
        for t, directory in search_paths:
            localdict = dict()
            #Get the files list
            file_list = glob(path.join(directory, self._glob_pattern))
            #Start from the run day
            for _ in range(10):
                for filename in file_list:
                    #Get the base name
                    bn = path.basename(filename)
                    datestring = t.strftime("%Y%m%d")
                    if is_phys:
                        if (bn.find(t.strftime("%Y%m%d_f_")) != -1) and self._time_interval.contains(t):
                            #Build the tuple and append the tuple to the outputs
                            if localdict.has_key(datestring):
                                localdict[datestring].append(filename)
                            else:
                                localdict[datestring]=[filename]
                    else:
                        if (bn.find(t.strftime("%Y%m%d")) != -1) and self._time_interval.contains(t):
                            localdict[datestring]=[filename]
                t += timedelta(1)
            wdict[directory] = localdict

        UNIQUE_dict=dict()
        for t, directory in search_paths:
            localdict=wdict[directory]
            for key in localdict.keys():
                UNIQUE_dict[key] = localdict[key] # this can overwrite

        output = list()
        key_list = sorted(UNIQUE_dict.keys())
        for t in key_list:
            filelist=UNIQUE_dict[t]
            for filename in filelist:
                output.append((t,filename))

        return output



    def extract_analysis(self, outputdir, rundays=[2], command="gzip -cd $INFILE > $OUTFILE", remove_ext=True):
        """
        Extracts analysis files to outputdir.

        Args:
            - *outputdir*: path to the output directory. If it does not exists
              this method will attempt to create it. If it exists but it is not
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
        outputdir = path.abspath(str(outputdir))
        return self._extract(file_list, outputdir, str(command), remove_ext)

    def extract_forecast(self, outputdir, rundays=[2,5], command="gzip -cd $INFILE > $OUTFILE", remove_ext=True):
        """
        Extracts forecast files to outputdir.

        Args: see extract_analysis.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        file_list = self.get_forecast_days(rundays)
        outputdir = path.abspath(str(outputdir))
        return self._extract(file_list, outputdir, str(command), remove_ext)
    
    def extract_simulation(self, outputdir, rundays=[2], command="gzip -cd $INFILE > $OUTFILE", remove_ext=True):
        """
        Extracts forecast files to outputdir.

        Args: see extract_analysis.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        file_list = self.get_simulation_days(rundays)
        outputdir = path.abspath(str(outputdir))
        return self._extract(file_list, outputdir, str(command), remove_ext)

    def extract_from_list(self, file_list, outputdir, command="gzip -cd $INFILE > $OUTFILE", remove_ext=True):
        """
        Extract the files from a list into outputdir.

        Args:
            - *file_list*: a list of tuples (datetime, filename).
            - *outputdir*: path to the output directory. If it does not exists
              this method will attempt to create it. If it exists but it is not
              a directory a ValueError will be raised.
            - *command* (optional): the command string template. See extract_analysis.
            - *remove_ext* (optional): boolean value. See extract_analysis.

        Returns: a list of paths pointing to the successfully extracted files.
        """
        return self._extract(file_list, outputdir, str(command), remove_ext)

    #Private methods
    def _extract(self, file_list, outputdir, command, remove_ext):
        #If outputdir exists make sure it's a directory
        if path.exists(outputdir):
            if not path.isdir(outputdir):
                raise ValueError("%s is not a directory" % outputdir)
        else:
            #Try to create the output directory
            mkdir(outputdir)
        #Build command template
        if (command.find('$INFILE') != -1) and (command.find('$OUTFILE') != -1):
            template = Template(command)
        else:
            raise ValueError("Invalid command string: " + command)
        output = list()
        #For each input file
        for _,ifn in file_list:
            #Get the path for the output file
            ofn = path.basename(ifn)
            if remove_ext:
                ofn = path.splitext(ofn)[0]
            ofn = path.join(outputdir, ofn)
            #Build the actual command string
            comstring = template.substitute(INFILE=ifn, OUTFILE=ofn)
            #Lauch the command
            retcode = SP.call(comstring, shell=True)
            #If the subprocess ended normally
            if retcode == 0:
                output.append(ofn)
            else:
                warn("Extract command exit code: %d" % retcode)
        return output

    #Static Methods
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
        sat_archive = str(sat_archive)
        if path.exists(sat_archive):
            if not path.isdir(sat_archive):
                raise ValueError("%s is not a directory" % sat_archive)
        else:
            raise ValueError("%s doesn't exist" % sat_archive)
        sat_files = glob(path.join(sat_archive, glob_pattern))
        output = list()
        for dt,p in L:
            ds = dt.strftime("%Y%m%d")
            for f in sat_files:
                bn = path.basename(f)
                if bn.find(ds) != -1:
                    output.append((dt,p,f))
        return output
