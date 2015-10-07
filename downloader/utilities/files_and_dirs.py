# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from os import makedirs
from os.path import exists, isdir

class DirIsAFileError(IOError):
    """An error that means that a path to a dir points to a file"""
    pass

def ensure_dir(path, log, expected=False):
    """
    Check if a directory exists and, if not, create it. If needed, create
    also all the intermediate directory between the last existing folder
    and the requested one.
    
    Args:
        - *path*: The path of the desidered directory as a string
        - *log*: an object of the Log class to report the operations done
          during the execution
        - *expected* (optional): If expected is True then, if the directory
          is missing, it will be reported to the user in a severity level 
          that is not the debug one but it is higher. Default is True
    
    Raises:
        - *DirIsAFileError* if the path points to a file    
    """
    if exists(path):
        if isdir(path):
           log.debug(path + " exists and it is a directory!")
        else:
            log.error(path + " exists but it is not a"
                      "directory!")
            raise DirIsAFileError
    else:
        if expected:
            log.info("Dir " + path + " not found. Now it will "
                     "be created!")
        makedirs(path)
        log.debug("Directory " + path + " created")

