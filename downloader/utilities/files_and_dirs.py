from os import makedirs
from os.path import exists, isdir

class DirIsAFileError(IOError):
    pass

def ensure_dir(path, log, expected=False):
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

