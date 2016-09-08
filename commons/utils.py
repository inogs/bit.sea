# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from __future__ import print_function
import os,sys
import numpy as np
import re

"""Helper functions"""

def is_number(val):
    """Tells if a value is a number type

    Args:
        - *val*: the value to test

    Returns:
        - True if val is a number type.
        - False if val is not a number type.
    """
    return isinstance(val, (int, long, float, complex, np.float32))

def get_date_string(s):
    """Finds a date in YYYYMMDD format into a bigger string.

    Args:
        - *s*: the string.

    Returns: a tuple with the date string in YYYY-MM-DD format as first element
    and YYYYMMDD format as second element. If no date has been found returns
    two empty strings.
    """
    m = re.search('.*([0-9]{4})([0-9]{2})([0-9]{2}).*',s)
    if m is None:
        return "", ""
    longdate = "%s-%s-%s" % (m.group(1), m.group(2), m.group(3))
    shortdate = "%s%s%s" % (m.group(1), m.group(2), m.group(3))
    return longdate, shortdate

def addsep(string):
    if string[-1] != os.sep:
        return string + os.sep
    else:
        return  string

def file2stringlist(filename):
    '''
    Argument : string - a path name indicating a text file
    Returns  : a list of strings
  
    A text file is converted in a list of strings, one for each row.
    
    '''
    LIST=[]
    filein=open(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST

def find_index(thestring, STRINGLIST):
    '''
    Searches a string in a list of strings.
    The string is supposed exist in the list. If not, the function raises an error.
    The string list is supposed be like in NetCDF files, where a list of variable names is
    written in an only 2D character array. 
    '''
    nStrings = STRINGLIST.shape[0]
    for istring in range(nStrings):
        strippedstring=STRINGLIST[istring,:].tostring().strip()
        if strippedstring == thestring: break
    else:
        print(thestring + " Not Found")
        raise NameError('Variable should be one of the following: ' + str([STRINGLIST[istring,:].tostring().strip() for istring in range(nStrings)]))
    return istring

def die(why, exit_code=1, print_usage=True):
    print("FATAL ERROR: " +  str(why), file=sys.stderr)
    sys.exit(exit_code)

def is_valid_path(path, is_dir_check=False):
    if os.path.exists(path):
        if is_dir_check:
            if os.path.isdir(path):
                return path
            else:
                die("'%s' is not a directory." % (path,))
        else:
            return path
    else:
        die("'%s' is not a valid path." % (path,))
