# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def mapplot(map_dict, fig, ax, mask=None, min_ticks=4, max_ticks=8, cbar_ticks=5, coastline=False, dpi=72.0):
    """Map plotting procedure (draft)

    Args:
        - *map_dict*: a dictionary as built by get_maps_data method of
          MapBuilder.
        - *mask* (optional): a Mask object that will be used to set the ticks.
        - *min_ticks* (optional): Number of ticks to set in the shorter axis.
        - *max_ticks* (optional): Number of ticks to set in the longer axis.
    """
    shape = map_dict['data'].shape
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
    ax.clear()
    fig.set_dpi(dpi)
    fig.set_size_inches(shape[1] / dpi, shape[0] / dpi)
    if (not mask is None) and coastline:
        coast_m = np.array(mask.mask[0,:,:], dtype=np.float32)
        coast_m[coast_m != 0] = np.nan
        ax.imshow(coast_m, cmap='hot')
        ax.invert_yaxis()
        fig.hold(True)
    clim = map_dict['clim']
    title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
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
        x_points = np.linspace(0,shape[1]-1,x_ticks).tolist()
        x_labels = list()
        for x in x_points:
            x_labels.append(int(round(mask.convert_x_y_to_lon_lat(int(x), 0)[0])))
        ax.set_xticks(x_points)
        ax.set_xticklabels(x_labels)
        #Set Y axis ticks
        y_points = np.linspace(0,shape[0]-1,y_ticks).tolist()
        y_labels = list()
        for y in y_points:
            y_labels.append(int(round(mask.convert_x_y_to_lon_lat(0, int(y))[1])))
        ax.set_yticks(y_points)
        ax.set_yticklabels(y_labels)
    fig.suptitle(title)
    ax.grid()
    return fig, ax
