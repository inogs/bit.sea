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
    parser.add_argument(   '--logo', '-l',
                                type = str,
                                required = False,
                                default = "/galileo/home/userexternal/gbolzon0/LogoEchoOGS4.png",
                                help = 'logo file')
    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')
import sys
from commons.utils import is_valid_path, nan_compare
import numpy as np
from glob import glob
from commons.Timelist import TimeInterval, TimeList
from commons.mask import Mask
from mpl_toolkits.basemap import Basemap
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

xlim=[-6.5,36.5]
ylim=[30,46]
xC=(xlim[0]+xlim[1])/2
yC=(ylim[0]+ylim[1])/2
 
mapobj = Basemap(projection='merc',lat_0=xC,lon_0=yC,\
                                  llcrnrlon = xlim[0], \
                                  llcrnrlat = ylim[0], \
                                  urcrnrlon = xlim[1], \
                                  urcrnrlat = ylim[1], \
                                  area_thresh=None, \
                                  resolution='i')



TI = TimeInterval("1950","2050","%Y")
TL = TimeList.fromfilenames(TI, inputdir, file_pattern )
TheMask = Mask(maskfile)
mb = MapBuilder(plotlistfile, TL, TheMask, outputdir)

logo=mb.read_background(args.logo)
mb.plot_maps_data(mapobj,background_img=logo, maptype=1, nranks=nranks, rank=rank)

