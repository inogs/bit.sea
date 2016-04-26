import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Plots profile images from NetCDF  files
    The schema is : Float vs Forecast, Analysis
    '''
    ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--forecast_dir','-f',
                                type = str,
                                required = True,
                                help = 'directory of forecast NetCDF files')
    parser.add_argument(   '--analysis_dir','-a',
                                type = str,
                                required = True,
                                help = 'directory of forecast NetCDF files')

    return parser.parse_args()

args = argument()

import pylab as pl
import numpy as np
import matplotlib.patches as mpatches
import scipy.io.netcdf as NC
import libxmp, libxmp.utils
from libxmp import XMPFiles, consts
from layer_integral import coastline

def figure_generator(p):
    ''' Generates a figure to plot the matchups related to a bioFloat cycle
    There are 6 axes: map, temperature, salinity, chl, oxygen and nitrate

    Arguments:
    * p * is a profile object

    Returns
    fig, axes (array of axes handlers)
    '''
    fig, axs = pl.subplots(2,3, facecolor='w', edgecolor='k')
    hsize=10
    vsize=12
    fig.set_size_inches(hsize,vsize)
    #figtitle = " date="+p.time.strftime('%Y/%m/%d')+" float="+p.name()
    #fig.set_title(figtitle)
    fig.subplots_adjust(hspace = 0.15, wspace=0.3)
    axs = axs.ravel()

    ax = axs[0]
    c_lon, c_lat=coastline.get()
    ax.plot(c_lon,c_lat, color='#000000',linewidth=0.5)
    ax.plot(p.lon,p.lat,'ro')
    ax.set_xticks(np.arange(-6,36,2))
    ax.set_yticks(np.arange(0,100,2))
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(p.time.strftime('%Y/%m/%d'))
    extent=10 #degrees
    ax.set_xlim([p.lon -extent/2, p.lon+extent/2])
    ax.set_ylim([p.lat -extent/2, p.lat+extent/2])
    bbox=ax.get_position()

    deltax, _ =bbox.size
    new_deltay = deltax* hsize/vsize
    bottom = bbox.ymax - new_deltay
    ax.set_position([bbox.xmin, bottom, deltax, new_deltay])

    floatlabel = 'Float \n'+ p.name() +" - "+str(p._my_float.cycle)
    b_patch = mpatches.Patch(color='red', label='Model')
    g_patch = mpatches.Patch(color='blue', label=floatlabel)
    ax.legend(handles=[b_patch,g_patch], bbox_to_anchor=(0, -0.5), loc=2)

    for ax in axs[1:]:
        ax.set_ylim(0,400)
        ax.locator_params(axis='x',nbins=4)
        ax.yaxis.grid()

    for ax in [axs[2], axs[4], axs[5]]:
        ax.set_yticklabels([])

    return fig,axs



def ncreader(filenc,zlevels_out,profileobj):
    ''' 
    '''

    depths = len(zlevels_out)
    f = NC.netcdf_file(filenc, 'r')
    f.close()



def add_metadata(filepng,p):
    xmpfile = XMPFiles( file_path=filepng, open_forupdate=True )
    xmp = xmpfile.get_xmp()
    if xmp is None:
        xmp = libxmp.XMPMeta()
    xmp.set_property(consts.XMP_NS_DC, 'float', p.name() )
    xmp.set_property(consts.XMP_NS_DC, 'date', p.time.strftime('%Y%m%d') )
    xmp.set_property(consts.XMP_NS_DC, 'hour', p.time.strftime('%H:%M:%S') )
    xmp.set_property(consts.XMP_NS_DC, 'position.lat',str(p.lat)+"N")
    xmp.set_property(consts.XMP_NS_DC, 'position.lon',str(p.lon)+"E")
    xmpfile.put_xmp(xmp)
    xmpfile.close_file()
