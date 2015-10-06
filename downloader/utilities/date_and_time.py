# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from time import strftime

def now_as_string():
    '''Return the current date as a string in the format "YYYYMMDD"'''
    return strftime("%Y%m%d")

