#!/usr/bin/env python
# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
# Script to build layer maps from model files

from __future__ import print_function
import sys
import os
from glob import glob

try:
    from layer_integral.mapbuilder import MapBuilder
except ImportError:
    print("You should run this script from the bit.sea root directory.", file=sys.stderr)
    sys.exit(2)

###Default parameters###
#Input directory
inputdir = "/pico/home/usera07ogs/a07ogs00/OPA/V4/wrkdir/2/POSTPROC/AVE_FREQ_1/"

#Output directory
outputdir = "/pico/scratch/usera07ogs/a07ogs02/layer_maps"

#File pattern for glob
file_pattern = "ave*nc"

#Path to plot list XML file
plotlistfile = "/pico/home/usera07ogs/a07ogs02/bit.sea/postproc/Plotlist.xml"

#Path to mask file
maskfile = "/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc"

def usage():
    print("""
Usage:
    %s [inputdir] [outputdir] [plotlist] [maskfile]

""" % (sys.argv[0]))
    print("Default input directory: '%s'" % (str(inputdir),))
    print("Default output directory: '%s'" % (str(outputdir),))
    print("Default plot list file: '%s'" % (str(plotlistfile),))
    print("Default mask file: '%s'" % (str(maskfile),))

def die(why, exit_code=1):
    print("FATAL ERROR: " +  str(why), file=sys.stderr)
    usage()
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

def validate_arguments(arg1,arg2,arg3,arg4):
    inputdir = is_valid_path(arg1, True)
    outputdir = is_valid_path(arg2, True)
    plotlistfile = is_valid_path(arg3)
    maskfile = is_valid_path(arg4)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        #No arguments, let's validate the defaults
        pass
    elif len(sys.argv) == 2:
        #Got one argument
        inputdir = sys.argv[1]
        validate_arguments(inputdir, outputdir, plotlistfile, maskfile)
    elif len(sys.argv) == 3:
        #Got two arguments
        inputdir = sys.argv[1]
        outputdir = sys.argv[2]
    elif len(sys.argv) == 4:
        #Got three arguments
        inputdir = sys.argv[1]
        outputdir = sys.argv[2]
        plotlistfile = sys.argv[3]
    elif len(sys.argv) == 5:
        #Got four arguments
        inputdir = sys.argv[1]
        outputdir = sys.argv[2]
        plotlistfile = sys.argv[3]
        maskfile = sys.argv[4]
    else:
        usage()
        sys.exit(1)

    try:
        validate_arguments(inputdir, outputdir, plotlistfile, maskfile)
        file_list = glob(inputdir + "/" + file_pattern)
        mb = MapBuilder(plotlistfile, file_list, maskfile, outputdir)
        mb.plot_maps_data()
    except Exception as e:
        die(e)
