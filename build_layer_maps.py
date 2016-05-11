#!/usr/bin/env python
# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>


from __future__ import print_function
import sys
import os
import numpy as np
from glob import glob
import argparse
import time
from layer_integral import coastline

def argument():
    parser = argparse.ArgumentParser(description = '''
    Script to build layer maps from model files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =False,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V4/wrkdir/2/POSTPROC/AVE_FREQ_1/TMP/",
                                help = ''' Input directory'''

                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                default = '/pico/scratch/usera07ogs/a07ogs02/layer_maps/',
                                help = 'Output directory')

    parser.add_argument(   '--pattern', '-g',
                                type = str,
                                required = False,
                                default = 'ave*nc',
                                help = 'glob search pattern')

    parser.add_argument(   '--plotlistfile', '-p',
                                type = str,
                                required = False,
                                default = 'postproc/Plotlist.xml',
                                help = 'Path to plot list XML file')

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = False,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                help = 'Path to mask file')
    return parser.parse_args()

try:
    from layer_integral.mapbuilder import MapBuilder
except ImportError:
    print("You should run this script from the bit.sea root directory.", file=sys.stderr)
    sys.exit(2)

args = argument()




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


maskfile     = is_valid_path(args.maskfile)
plotlistfile = is_valid_path(args.plotlistfile)
inputdir     = is_valid_path(args.inputdir,True)
outputdir    = is_valid_path(args.outputdir,True)
file_pattern = args.pattern



try:
    c_lon, c_lat=coastline.get()
    # Elimination of some parts of the coastline,
    # in order to leave more space for text, if needed
    ii = (c_lat > 40.0) & (c_lon < 0.0) # atlantic coast
    jj = (c_lat > 42.0) & (c_lon > 26 ) # black sea
    c_lon[ii | jj] = np.NaN
    c_lat[ii | jj] = np.NaN
    ii = (c_lon < -6) | (c_lon > 36 )# out of box
    c_lon[ii] = np.NaN
    c_lat[ii] = np.NaN

except:
    c_lon=None
    c_lat=None

try:
    file_list = glob(inputdir + "/" + file_pattern)
    mb = MapBuilder(plotlistfile, file_list, maskfile, outputdir)
    mb.plot_maps_data(coastline_lon=c_lon, coastline_lat=c_lat)
    mb.plot_maps_data(coastline_lon=c_lon, coastline_lat=c_lat,maptype=1)
    mb.plot_maps_data(coastline_lon=c_lon, coastline_lat=c_lat,maptype=2)
except Exception as e:
    die(e, 2, False)


