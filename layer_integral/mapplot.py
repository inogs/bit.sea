# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import code
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def mapplot(map_dict, fig, ax, mask=None,ncolors=256,cbar_ticks=5, coastline_lon=None, coastline_lat=None, dpi=72.0):
    """Map plotting procedure (draft)
    Hardcoded features:
        - colormap jet
        - watermark
        - ticks
        - size of the figure

    Args:
        - *map_dict*: a dictionary as built by get_maps_data method of
          MapBuilder.
        - *fig*: a reference to a Figure object, if None mapplot will create a new Figure.
        - *ax*: a reference to an Axes object, if None mapplot will create a new Figure.
        - *mask* (optional): a Mask object that will be used to set the ticks.
        - *ncolors* (optional) : the number of colors of colormap
        - *cbar_ticks* (optional): Number of ticks on the colorbar (default: 5).
        - *coastline_lon* (optional): Numpy array defining the coastline longitudes.
        - *coastline_lat* (optional): Numpy array defining the coastline latitudes.
        - *dpi* (optional): sets the DPI (default: 72.0).
    Returns:
        A figure and an Axes object that can be passed again to mapplot
    """
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
        fig.set_dpi(dpi)
        #shape = map_dict['data'].shape
        #fig.set_size_inches(shape[1] / float(dpi), shape[0] / float(dpi))
        fig.set_size_inches(10.0, 10.0*16/42)
    else:
        fig.clf()
        fig.add_axes(ax)
    ax.set_position([0.07, 0.11, 0.78, 0.85])
    clim = map_dict['clim']

    if not(mask is None):
        lon_min = mask.xlevels.min()
        lon_max = mask.xlevels.max()
        lat_min = mask.ylevels.min()
        lat_max = mask.ylevels.max()
        cmap=pl.get_cmap('jet',ncolors)
        im = ax.imshow(map_dict['data'], extent=[lon_min, lon_max, lat_max, lat_min], cmap=cmap)
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
        x_points = np.arange(-6,36,4).tolist()
        y_points = np.arange(32,46,4).tolist()
        #Set X axis ticks
        ax.set_xticks(x_points)
        #Set Y axis ticks
        ax.set_yticks(y_points)

        if not ((coastline_lon is None) or (coastline_lat is None)):
            # Flatten coastline arrays
            coastline_lon = np.ravel(coastline_lon)
            coastline_lat = np.ravel(coastline_lat)
            if len(coastline_lon) != len(coastline_lat):
                raise ValueError("coastline arrays must have the same length")
            #Draw coastline
            ax.plot(coastline_lon,coastline_lat, color='#000000',linewidth=0.5)
            ax.set_xlim([lon_min, lon_max])
            ax.set_ylim([lat_min, lat_max])
    #ax.text(-7,44,map_dict['layer'].__repr__()  ,ha='left',va='center')
    #ax.text(-7,42,map_dict['date']   ,ha='left',va='center')
    #ax.text(-7,40,map_dict['varname'],ha='left',va='center')
    ax.set_xlabel('longitude (deg)')
    ax.set_ylabel('latitude (deg)')

    # watermark
    #ax.text(35, 46, 'OGS Echo Group',
    #     fontsize=60, color='gray',
    #     ha='right', va='top', alpha=0.3, rotation=18)

    title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
    fig.suptitle(title)

    return fig, ax


def mapplot_onlycolor(map_dict, fig, ax, mask=None,ncolors=256,cbar_ticks=5, coastline_lon=None, coastline_lat=None, dpi=72.0):
    """Map plotting procedure only map color(draft)
    Hardcoded features:
        - colormap jet
        - watermark
        - ticks
        - size of the figure

    Args:
        - *map_dict*: a dictionary as built by get_maps_data method of
          MapBuilder.
        - *fig*: a reference to a Figure object, if None mapplot will create a new Figure.
        - *ax*: a reference to an Axes object, if None mapplot will create a new Figure.
        - *mask* (optional): a Mask object that will be used to set the ticks.
        - *ncolors* (optional) : the number of colors of colormap
        - *cbar_ticks* (optional): Number of ticks on the colorbar (default: 5).
        - *coastline_lon* (optional): Numpy array defining the coastline longitudes.
        - *coastline_lat* (optional): Numpy array defining the coastline latitudes.
        - *dpi* (optional): sets the DPI (default: 72.0).
    Returns:
        A figure and an Axes object that can be passed again to mapplot
    """

    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
        fig.set_dpi(dpi)
        #shape = map_dict['data'].shape
        #fig.set_size_inches(shape[1] / float(dpi), shape[0] / float(dpi))
        fig.set_size_inches(10.0, 10.0*16/42)
    else:
        fig.clf()
        fig.add_axes(ax)
    ax.set_position([0.07, 0.11, 0.78, 0.85])
    clim = map_dict['clim']

    if not(mask is None):
        lon_min = mask.xlevels.min()
        lon_max = mask.xlevels.max()
        lat_min = mask.ylevels.min()
        lat_max = mask.ylevels.max()
        cmap=pl.get_cmap('jet',ncolors)
        im = ax.imshow(map_dict['data'], extent=[lon_min, lon_max, lat_max, lat_min], cmap=cmap)
    else:
        im = ax.imshow(map_dict['data'])
    #Set color bar
    im.set_clim(clim[0], clim[1])
    # cbar_ticks_list = np.linspace(clim[0], clim[1], cbar_ticks).tolist()
    # cbar_ticks_labels = list()
    # for t in cbar_ticks_list:
    #     cbar_ticks_labels.append("%g" % (t,))
    # div = make_axes_locatable(ax)
    # cax = div.append_axes("right", size="3%", pad=0.05)
    # cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks_list)
    # cbar.ax.set_yticklabels(cbar_ticks_labels)
    ax.invert_yaxis()
    if not mask is None:
        x_points = np.arange(-6,36,4).tolist()
        y_points = np.arange(32,46,4).tolist()
        #Set X axis ticks
        #ax.set_xticks(x_points)
        #Set Y axis ticks
        #ax.set_yticks(y_points)

        if not ((coastline_lon is None) or (coastline_lat is None)):
            # Flatten coastline arrays
            coastline_lon = np.ravel(coastline_lon)
            coastline_lat = np.ravel(coastline_lat)
            if len(coastline_lon) != len(coastline_lat):
                raise ValueError("coastline arrays must have the same length")
            #Draw coastline
            #ax.plot(coastline_lon,coastline_lat, color='#000000',linewidth=0.5)
            ax.set_xlim([lon_min, lon_max])
            ax.set_ylim([lat_min, lat_max])
    #ax.text(-7,44,map_dict['layer'].__repr__()  ,ha='left',va='center')
    #ax.text(-7,42,map_dict['date']   ,ha='left',va='center')
    #ax.text(-7,40,map_dict['varname'],ha='left',va='center')
    #ax.set_xlabel('longitude (deg)')
    #ax.set_ylabel('latitude (deg)')


    title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
    #fig.suptitle(title)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_axis_off()

    #pl.axis('off')
    #pl.show()
    return fig, ax


################################################################################

def mapplot_nocolor(map_dict, fig, ax, mask=None,ncolors=256,cbar_ticks=5, coastline_lon=None, coastline_lat=None, dpi=72.0):
    """Map plotting procedure only map(draft)
    Hardcoded features:
        - colormap jet
        - watermark
        - ticks
        - size of the figure

    Args:
        - *map_dict*: a dictionary as built by get_maps_data method of
          MapBuilder.
        - *fig*: a reference to a Figure object, if None mapplot will create a new Figure.
        - *ax*: a reference to an Axes object, if None mapplot will create a new Figure.
        - *mask* (optional): a Mask object that will be used to set the ticks.
        - *ncolors* (optional) : the number of colors of colormap
        - *cbar_ticks* (optional): Number of ticks on the colorbar (default: 5).
        - *coastline_lon* (optional): Numpy array defining the coastline longitudes.
        - *coastline_lat* (optional): Numpy array defining the coastline latitudes.
        - *dpi* (optional): sets the DPI (default: 72.0).
    Returns:
        A figure and an Axes object that can be passed again to mapplot
    """
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
        fig.set_dpi(dpi)
        fig.set_size_inches(10.0, 10.0*16/42)
    else:
        fig.clf()
        fig.add_axes(ax)
    ax.set_position([0.07, 0.11, 0.78, 0.85])
    clim = map_dict['clim']

    if not(mask is None):
        lon_min = mask.xlevels.min()
        lon_max = mask.xlevels.max()
        lat_min = mask.ylevels.min()
        lat_max = mask.ylevels.max()
        cmap=pl.get_cmap('jet',ncolors)
        im = ax.imshow(map_dict['data'], extent=[lon_min, lon_max, lat_max, lat_min], cmap=cmap)
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
    cbar.ax.set_yticklabels([])
    #cbar.ax.get_yticklabels().set_visible(False)
    ax.invert_yaxis()
    if not mask is None:
        x_points = np.arange(-6,36,4).tolist()
        y_points = np.arange(32,46,4).tolist()
        #Set X axis ticks
        ax.set_xticks(x_points)
        #Set Y axis ticks
        ax.set_yticks(y_points)

        if not ((coastline_lon is None) or (coastline_lat is None)):
            # Flatten coastline arrays
            coastline_lon = np.ravel(coastline_lon)
            coastline_lat = np.ravel(coastline_lat)
            if len(coastline_lon) != len(coastline_lat):
                raise ValueError("coastline arrays must have the same length")
            #Draw coastline
            ax.plot(coastline_lon,coastline_lat, color='#000000',linewidth=0.5)
            ax.set_xlim([lon_min, lon_max])
            ax.set_ylim([lat_min, lat_max])
    #ax.text(-7,44,map_dict['layer'].__repr__()  ,ha='left',va='center')
    #ax.text(-7,42,map_dict['date']   ,ha='left',va='center')
    #ax.text(-7,40,map_dict['varname'],ha='left',va='center')
    #ax.set_xlabel('longitude (deg)')
    #ax.set_ylabel('latitude (deg)')

    # watermark
    #ax.text(35, 46, 'OGS Echo Group',
    #     fontsize=60, color='gray',
    #     ha='right', va='top', alpha=0.3, rotation=18)

    title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
    #fig.suptitle(title)
    im.set_visible(False)
    return fig, ax
