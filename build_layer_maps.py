#!/usr/bin/env python
# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

from __future__ import print_function
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Script to build layer maps from model files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =False,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/wrkdir/2/POSTPROC/AVE_FREQ_1/TMP/",
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
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                help = 'Path to mask file')
    parser.add_argument(   '--background', '-b',
                                type = str,
                                required = False,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/etc/static-data/POSTPROC/background_medeaf.png",
                                help = 'Path to mask file')
    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')
import sys
from commons.utils import is_valid_path, nan_compare
import numpy as np
from glob import glob
from layer_integral import coastline
from commons.Timelist import TimeInterval, TimeList
from commons.mask import Mask
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False



try:
    from layer_integral.mapbuilder import MapBuilder
except ImportError:
    print("You should run this script from the bit.sea root directory.", file=sys.stderr)
    sys.exit(2)


def die(why, exit_code=1, print_usage=True):
    print("FATAL ERROR: " +  str(why), file=sys.stderr)
    sys.exit(exit_code)


maskfile     = is_valid_path(args.maskfile)
plotlistfile = is_valid_path(args.plotlistfile)
inputdir     = is_valid_path(args.inputdir,True)
outputdir    = is_valid_path(args.outputdir,True)
file_pattern = args.pattern

try:
    c_lon, c_lat=coastline.get()
    # Elimination of some parts of the coastline,
    # in order to leave more space for text, if needed
    ii = nan_compare(c_lat, ">", 40.0) & nan_compare(c_lon,"<", 0.0) # atlantic coast
    jj = nan_compare(c_lat, ">", 42.0) & nan_compare(c_lon, ">", 26 ) # black sea
    c_lon[ii | jj] = np.NaN
    c_lat[ii | jj] = np.NaN
    ii = nan_compare(c_lon, "<", -6) | nan_compare(c_lon, ">", 36)# out of box
    c_lon[ii] = np.NaN
    c_lat[ii] = np.NaN

except:
    c_lon=None
    c_lat=None

TI = TimeInterval("1950","2050","%Y")
TL = TimeList.fromfilenames(TI, inputdir, file_pattern )
TheMask = Mask(maskfile)
mb = MapBuilder(plotlistfile, TL, TheMask, outputdir)
#mb.plot_maps_data(coastline_lon=c_lon, coastline_lat=c_lat)
background=mb.read_background(args.background)
mb.plot_maps_data(coastline_lon=c_lon, coastline_lat=c_lat,background_img=background, maptype=1, nranks=nranks, rank=rank)
#mb.plot_maps_data(coastline_lon=c_lon, coastline_lat=c_lat,maptype=2)

