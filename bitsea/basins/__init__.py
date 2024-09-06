# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from basins.basin import Basin

from importlib import import_module
from inspect import getmembers

class BasinNotFoundError(Exception): pass

class get_basin(object):
    def __init__(self):
        self.__default = "OGS"

    def __call__(self, name, module=None):
        if module is None:
            module = self.__default
        else:
            self.__default = module

        mod = import_module('basins.' + module)
        for obj_name, obj in getmembers(mod):
            if isinstance(obj, Basin):
                if obj.name == name:
                    return obj
                elif obj.extended_name == name:
                    return obj
                else:
                    pass
        raise BasinNotFoundError 

    def set_default(self, default):
        self.__default = default

get = get_basin()
