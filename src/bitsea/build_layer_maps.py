#!/usr/bin/env python
# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import argparse
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import existing_file_path

def argument():
    parser = argparse.ArgumentParser(description = '''
    Script to build layer maps from model files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                required =False,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/wrkdir/2/POSTPROC/AVE_FREQ_1/TMP/",
                                help = ''' Input directory'''
                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                default = '/pico/scratch/usera07ogs/a07ogs02/layer_maps/',
                                help = 'Output directory')

    parser.add_argument(   '--pattern', '-g',
                                type = str,
                                required = False,
                                default = 'ave*nc',
                                help = 'glob search pattern')

    parser.add_argument(   '--plotlistfile', '-p',
                                type = existing_file_path,
                                required = False,
                                default = 'postproc/Plotlist.xml',
                                help = 'Path to plot list XML file')

    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
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
from bitsea.commons.utils import is_valid_path
from bitsea.commons.Timelist import TimeInterval, TimeList
from bitsea.commons.mask import Mask
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
    from bitsea.layer_integral.mapbuilder import MapBuilder
except ImportError:
    print("You should run this script from the bit.sea root directory.", file=sys.stderr)
    sys.exit(2)



maskfile     = args.maskfile
plotlistfile = args.plotlistfile
inputdir     = args.inputdir
outputdir    = args.outputdir
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
TheMask = Mask.from_file(maskfile)
mb = MapBuilder(plotlistfile, TL, TheMask, outputdir)

logo=mb.read_background(args.logo)
mb.plot_maps_data(mapobj,background_img=logo, maptype=1, nranks=nranks, rank=rank)

