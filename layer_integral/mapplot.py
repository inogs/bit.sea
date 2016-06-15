# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.font_manager as font_manager
from matplotlib.font_manager import FontProperties


def mapplot(map_dict, fig, ax, mask=None,ncolors=256,cbar_ticks=5, coastline_lon=None, coastline_lat=None, dpi=72.0):
    """Map plotting procedure (draft)
    Hardcoded features:
        - colormap jet
        - ticks
        - size of the figure

    Args:
        - *map_dict     *: a dictionary as built by get_maps_data method of MapBuilder.
        - *fig          *: a reference to a Figure object, if None mapplot will create a new Figure.
        - *ax           *: a reference to an Axes object, if None mapplot will create a new Figure.
        - *mask         * (optional): a Mask object that will be used to set the ticks.
        - *ncolors      * (optional) : the number of colors of colormap
        - *cbar_ticks   * (optional): Number of ticks on the colorbar (default: 5).
        - *coastline_lon* (optional): Numpy array defining the coastline longitudes.
        - *coastline_lat* (optional): Numpy array defining the coastline latitudes.
        - *dpi          * (optional): sets the DPI (default: 72.0).
    Returns:
        A figure and an Axes object that can be passed again to mapplot
        
    Examples:
        fig, ax = mapplot({'data':Map2d, 'clim':[0,1]})
        fig, ax = mapplot({'data':Map2d, 'clim':[0,1]}, fig, ax )
                
        from layer_integral import coastline
        clon,clat = coastline.get()
        fig, ax = mapplot({'data':Map2d, 'clim':[0,1]}, fig, ax, coastline_lon=clon, coastline_lat=clat)
        
        from commons.mask import Mask
        TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')
        fig, ax = mapplot({'data':Map2d, 'clim':[0,1]}, fig, ax, mask=TheMask, coastline_lon=clon, coastline_lat=clat)
        fig, ax = mapplot({'data':Map2d, 'clim':[0,1], 'date':longdate, 'layer':l}, fig, ax, mask=TheMask, coastline_lon=clon, coastline_lat=clat)
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
    ax.set_position([0.08, 0.11, 0.78, 0.78])
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
            ax.set_xlim([-6, 36])
            ax.set_ylim([30, 46])
    #ax.text(-7,44,map_dict['layer'].__repr__()  ,ha='left',va='center')
    #ax.text(-7,42,map_dict['date']   ,ha='left',va='center')
    #ax.text(-7,40,map_dict['varname'],ha='left',va='center')
    ax.set_position([0.08, 0.11, 0.78, 0.78])
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    #ax.set_xlabel('longitude (deg)')
    #ax.set_ylabel('latitude (deg)')
    # watermark
    #ax.text(35, 46, 'OGS Echo Group',
    #     fontsize=60, color='gray',
    #     ha='right', va='top', alpha=0.3, rotation=18)
    if map_dict.has_key('layer'):
        title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
    else:
        title = "%s %s" % (map_dict['date'], map_dict['varname'])
    #fig.suptitle(title)
    return fig, ax


def mapplot_medeaf(map_dict, fig, ax, mask=None,ncolors=256, background_img=None):
    """
    Designed for web site
    """
    
    #sfondo = pl.imread('/pico/home/userexternal/gbolzon0/bit.sea/layer_integral/20160610_OGS_slider_sito.png')

    font=FontProperties()
    font_prop   = font_manager.FontProperties(fname='TitilliumWeb-Regular.ttf', size=14)
    font_prop13 = font_manager.FontProperties(fname='TitilliumWeb-Regular.ttf', size=13)
    #font.set_name('TitilliumWeb')
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
        fig.set_size_inches(10.0, 10.0*16/42)
    else:
        fig.clf()
        fig.add_axes(ax)

    ratio=672./860
    ax.set_position([0.08, 0.13, ratio, ratio])
    clim = map_dict['clim']

    lon_min = mask.xlevels.min()
    lon_max = mask.xlevels.max()
    lat_min = mask.ylevels.min()
    lat_max = mask.ylevels.max()
    cmap=pl.get_cmap('jet',ncolors)
    im = ax.imshow(map_dict['data'], extent=[lon_min, lon_max, lat_max, lat_min], cmap=cmap)
    if not (background_img is None) : 
        ax.imshow(background_img, extent=[-6, 36, 30, 46])

    #Set color bar
    im.set_clim(clim[0], clim[1])
    cbar_ticks_list = np.linspace(clim[0], clim[1], 5).tolist()
    cbar_ticks_labels = list()
    for t in cbar_ticks_list:
        cbar_ticks_labels.append("%g" % (t,))


    cax = fig.add_axes((0.88,.13, 0.03, 0.78))
    cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks_list)
    cbar.ax.set_yticklabels(cbar_ticks_labels, 'fontproperties', font_prop)
    ax.invert_yaxis()

    ax.set_yticklabels([], 'fontproperties', font_prop)
    ax.set_xticklabels([], 'fontproperties', font_prop)
    ax.set_xticks(np.arange(-2,36,4).tolist())
    ax.set_yticks(np.arange(32,46,4).tolist())

    #ax.set_yticklabels([u"32N", u"36N", u"40N", u"44N"], 'fontproperties', font_prop)
    ax.set_xlim([-6, 36])
    ax.set_ylim([30, 46])
    t=ax.set_xlabel("longitude (deg)");
    t.set_font_properties(font_prop13)
    t=ax.set_ylabel("latitude (deg)")
    t.set_font_properties(font_prop13)



    ax.set_position([0.08, 0.13, ratio, ratio])
    ax.axes.get_xaxis().set_visible(True)
    ax.axes.get_yaxis().set_visible(True)
    units = map_dict['units']

    title_1 = "%s, %s "  % (map_dict['longname'],units)
    title_2 = "%s" %  map_dict['layer'].__repr__()
    title_3 = "%s" % (map_dict['date']).strftime('%d - %m - %Y')
    t1=ax.text(-6,46.3, title_1, verticalalignment='bottom', horizontalalignment='left')
    t2=ax.text(15,46.3, title_2, verticalalignment='bottom', horizontalalignment='center')
    t3=ax.text(36,46.3, title_3, verticalalignment='bottom', horizontalalignment='right')
    watermarkstring='Copyright : \nOGS ECHO GROUP\nmedeaf.inogs.it'
    ax.text(17, 35,watermarkstring ,fontsize=8,fontweight='bold', color='gray', ha='left', va='top', alpha=0.3)

    t1.set_font_properties(font_prop)
    t2.set_font_properties(font_prop)
    t3.set_font_properties(font_prop)
    return fig, ax


def mapplot_onlycolor(map_dict, fig, ax, mask=None,ncolors=256,cbar_ticks=5, dpi=72.0):
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
        - *dpi* (optional): sets the DPI (default: 72.0).
    Returns:
        A figure and an Axes object that can be passed again to mapplot
    """
    #watermark = pl.imread('/pico/home/userexternal/gbolzon0/griglia_senza_tacche.png')
    #watermark = pl.imread('/pico/home/userexternal/gbolzon0/watermark_ogs.png')
    #watermark = pl.imread('/pico/home/userexternal/gbolzon0/ogs_wm_72dpi.png')
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
    ax.set_position([0, 0, 1, 1])
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
    #cbar.set_visible(False)

    #ax.imshow(watermark,extent=[-4,9,31,35])
    #ax.imshow(watermark,extent=[lon_min,lon_max,lat_min,lat_max])
    watermarkstring='Copyright : \nOGS ECHO GROUP\nmedeaf.inogs.it'
    ax.text(-3, 33,watermarkstring ,fontsize=8,fontweight='bold', color='k', ha='left', va='top') # alpha=0.3)
    ax.text(30, 40,watermarkstring ,fontsize=8,fontweight='bold', color='k', ha='left', va='top')
    ax.text(17, 35,watermarkstring ,fontsize=8,fontweight='bold', color='gray', ha='left', va='top', alpha=0.3)
    ax.invert_yaxis()

    ax.set_xlim([-6, 36])
    ax.set_ylim([30, 46])

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_axis_off()

#     fig.savefig('try.jpg',dpi=72,quality=75)
#
#     r=pl.imread('try.jpg')
#     fig,ax=pl.subplots()
#     fig.set_size_inches(10.0, 10.0*16/42)
#     ax.set_position([0, 0, 1, 1])
#     ax.imshow(r)
#
#     ax.axes.get_xaxis().set_visible(False)
#     ax.axes.get_yaxis().set_visible(False)
#     ax.set_axis_off()
#     ax2 = fig.add_axes((0.03, -.05, 0.35, 0.4))
#     ax2.imshow(watermark)
#     ax2.axes.get_xaxis().set_visible(False)
#     ax2.axes.get_yaxis().set_visible(False)
#     ax2.set_axis_off()
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
            ax.set_xlim([-6, 36])
            ax.set_ylim([30, 46])
    #ax.text(-7,44,map_dict['layer'].__repr__()  ,ha='left',va='center')
    #ax.text(-7,42,map_dict['date']   ,ha='left',va='center')
    #ax.text(-7,40,map_dict['varname'],ha='left',va='center')
    #ax.set_xlabel('longitude (deg)')
    #ax.set_ylabel('latitude (deg)')

    # watermark
    #ax.text(35, 46, 'OGS Echo Group',
    #     fontsize=60, color='gray',
    #     ha='right', va='top', alpha=0.3, rotation=18)

    #title = "%s %s %s" % (map_dict['date'], map_dict['varname'], map_dict['layer'].__repr__())
    #fig.suptitle(title)
    im.set_visible(False)
    return fig, ax

if __name__ == '__main__':
    from commons.mask import Mask
    from commons.dataextractor import DataExtractor
    from datetime import datetime
    maskfile='/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc'
    mask = Mask(maskfile)
    filename='/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/ave.20000116-12:00:00.nc'

    DE = DataExtractor(mask,filename,'N1p')
    k=0
    map2d=DE.values[k,:,:]
    map2d[~mask.mask[k,:,:]] = np.NaN
    from commons.layer import Layer
    map_dict ={'data':map2d, 'clim':[0,0.1],'date':datetime.strptime('20160116','%Y%m%d'),'varname':'N1p', 'layer':Layer(0,10),'longname':'Phosphate', 'units':'mmolP/m3'}
    #fig, ax = mapplot_onlycolor(map_dict, fig=None, ax=None, mask=mask,ncolors=24,cbar_ticks=5, dpi=72.0)
    #fig.savefig('prova.jpg',dpi=72,quality=75)
    from layer_integral import coastline
    clon,clat = coastline.get()
    fig, ax = mapplot_medeaf(map_dict, fig=None, ax=None, mask=mask, ncolors=24)
    #(map_dict, fig=None, ax=None, mask=mask, coastline_lon=clon, coastline_lat=clat)
    #fig, ax = mapplot(map_dict, fig=None, ax=None, mask=mask,ncolors=24,cbar_ticks=5, dpi=72.0)
    #fig.show()
    fig.savefig('prova.png',dpi=86)

