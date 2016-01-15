# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def mapplot(map_dict, fig, ax, mask=None, min_ticks=4, max_ticks=8, cbar_ticks=5, coastline_lon=None, coastline_lat=None, dpi=72.0):
    """Map plotting procedure (draft)

    Args:
        - *map_dict*: a dictionary as built by get_maps_data method of
          MapBuilder.
        - *fig*: a reference to a Figure object, if None mapplot will create a new Figure.
        - *ax*: a reference to an Axes object, if None mapplot will create a new Figure.
        - *mask* (optional): a Mask object that will be used to set the ticks.
        - *min_ticks* (optional): Number of ticks to set in the shorter axis.
        - *max_ticks* (optional): Number of ticks to set in the longer axis.
        - *cbar_ticks* (optional): Number of ticks on the colorbar (default: 5).
        - *coastline_lon* (optional): Numpy array defining the coastline longitudes.
        - *coastline_lat* (optional): Numpy array defining the coastline latitudes.
        - *dpi* (optional): sets the DPI (default: 72.0).
    Returns:
        A figure and an Axes object that can be passed again to mapplot
    """
    shape = map_dict['data'].shape
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
        fig.set_dpi(dpi)
        fig.set_size_inches(shape[1] / float(dpi), shape[0] / float(dpi))
    else:
        fig.clf()
        fig.add_axes(ax)
    clim = map_dict['clim']
    title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
    if not(mask is None):
        lon_min = min(mask.xlevels)
        lon_max = max(mask.xlevels)
        lat_min = min(mask.ylevels)
        lat_max = max(mask.ylevels)
        im = ax.imshow(map_dict['data'], extent=[lon_min, lon_max, lat_max, lat_min])
    else:
        im = ax.imshow(map_dict['data'])
    #Set color bar
    im.set_clim(clim[0], clim[1])
    cbar_ticks_list = np.linspace(clim[0], clim[1], cbar_ticks).tolist()
    cbar_ticks_labels = list()
    for t in cbar_ticks_list:
        cbar_ticks_labels.append("%g" % (t,))
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="3%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks_list)
    cbar.ax.set_yticklabels(cbar_ticks_labels)
    ax.invert_yaxis()
    if not mask is None:
        #Decide who gets the most ticks
        if shape[0] > shape[1]:
            x_ticks = min_ticks
            y_ticks = max_ticks
        elif shape[1] > shape[0]:
            x_ticks = max_ticks
            y_ticks = min_ticks
        else:
            x_ticks = max_ticks
            y_ticks = max_ticks
        #Set X axis ticks
        x_points = np.linspace(lon_min,lon_max,x_ticks).tolist()
        ax.set_xticks(x_points)
        #Set Y axis ticks
        y_points = np.linspace(lat_min,lat_max,y_ticks).tolist()
        ax.set_yticks(y_points)
        if not ((coastline_lon is None) or (coastline_lat is None)):
            # Flatten coastline arrays
            coastline_lon = np.ravel(coastline_lon)
            coastline_lat = np.ravel(coastline_lat)
            if len(coastline_lon) != len(coastline_lat):
                raise ValueError("coastline arrays must have the same length")
            #Draw coastline
            ax.plot(coastline_lon,coastline_lat, color='#000000')
    fig.suptitle(title)
    ax.grid()
    return fig, ax
