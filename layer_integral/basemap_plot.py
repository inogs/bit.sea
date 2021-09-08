from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import numpy as np
from basins.region import Rectangle
import matplotlib.font_manager as font_manager
from mapplot import is_in_boxes, set_font, set_invisible

xlim=[-6.5,36.5]
ylim=[30,46]
xC=(xlim[0]+xlim[1])/2
yC=(ylim[0]+ylim[1])/2


map_obj = Basemap(projection='merc',lat_0=xC,lon_0=yC,\
                                  llcrnrlon = xlim[0], \
                                  llcrnrlat = ylim[0], \
                                  urcrnrlon = xlim[1], \
                                  urcrnrlat = ylim[1], \
                                  area_thresh=None, \
                                  resolution='i')


UNIT_LABELS_DICT=  {'fluxCO2'       : 'mmol C m'  + u'\u207B' + u'\u00B2' +'d' + u'\u207B' + u'\u00B9' ,\
                    'fluxco2'       : 'mg C m'    + u'\u207B' + u'\u00B2' +'d' + u'\u207B' + u'\u00B9' ,\
                    'mmolAlk/m3'    : 'mmol Eq m' + u'\u207B' + u'\u00B3'                              ,\
                    'mmolAlk/Kg'    : 'mmol Eq Kg'+ u'\u207B' + u'\u00B9'                              ,\
                    'mmolC/m3'      : 'mmol C m'  + u'\u207B' + u'\u00B3'                              ,\
                    'mmolC/Kg'      : 'mmol C Kg '+ u'\u207B' + u'\u00B9'                              ,\
                    'mgchla/m3'     : 'mg chla m' + u'\u207B' + u'\u00B3'                              ,\
                    'mgC/m3'        : 'mg C m'    + u'\u207B' + u'\u00B3'                              ,\
                    'mgC/m3/d'      : 'mg C m'    + u'\u207B' + u'\u00B3' +'d' + u'\u207B' + u'\u00B9' ,\
                    'gC/m3/y'       : 'g C m'     + u'\u207B' + u'\u00B3' +'y' + u'\u207B' + u'\u00B9' ,\
                    'mmolP/m3'      : 'mmol P  /m'+ u'\u00B3'                                          ,\
                    'mmol/m3'       : 'mmol /m'+ u'\u00B3'                                          ,\
                    'mmolN/m3'      : 'mmol N  /m'+ u'\u00B3'                                          ,\
                    'mmolSi/m3'     : 'mmol Si m' + u'\u207B' + u'\u00B3'                              ,\
                    'mmolO2/m3'     : 'mmol O'    + u'\u2082 '+' m' + u'\u207B' + u'\u00B3'            ,\
                    'm'             : 'm'                                                              ,\
                    'CO2'           : 'CO' + u'\u2082'                                                 }


plot_conf=np.dtype([('title'          ,'S100')  , ('var' ,'S100')        ,\
                    ('unit_conversion',np.float), ('unit_label','S100')  ,\
                    ('operator'       ,'S100')                           ,\
                    ('vmin'           ,np.float), ('vmax',np.float)      ,\
                    ('depth'          ,'S100')                           ,\
                    ('lonS'           ,np.float), ('lonE',np.float)      ,\
                    ('latS'           ,np.float), ('latE',np.float)      ,\
                    ('start'          ,'S100')  , ('end'   ,'S100')      ,\
                    ('indataAve'      ,'S100')  , ('aggregation' ,'S100'),\
                    ('InputDir'       ,'S100')  , ('OutputDir'   ,'S100'),\
                    ('FileName'       ,'S100')                           ])    

def pl_obj_convert(var, clim, starttime, endtime, units, layer_str):
    '''
    Prepares object for do_plot
    '''
    
    my_obj = np.ones((1,),dtype=plot_conf)
    my_obj['lonS']=-5.3
    my_obj['lonE']=36
    my_obj['latS']=30
    my_obj['latE']=46
    
    my_obj['title'] = var
    my_obj['vmin'],my_obj['vmax'] = clim
    my_obj['start'] = starttime
    my_obj['end']   = endtime
    my_obj['unit_label'] = units
    my_obj['depth'] = layer_str
    
    return my_obj[0]

def do_plot(myplot,matrixToPlot,maskObj):
# maskObj for coast line
    Lon = maskObj.xlevels
    Lat = maskObj.ylevels
    tk_m     = maskObj.getDepthIndex(0.5)
    tmask = maskObj.mask
    mask_2D   = tmask[tk_m,:,:]
    mask_2D[Lon < -5.3]   = False # Atlantic buffer
    Lonp=Lon#; Lonp[~mask_2D]=np.NaN
    Latp=Lat#; Latp[~mask_2D]=np.NaN

# Colors
    cmap=pl.get_cmap('jet')
    cmap=pl.get_cmap('gist_rainbow_r')
    cmap.set_bad(color='w',alpha=1.)
#-------------------------------------------------
    theLonMin=myplot['lonS']; theLonMax=myplot['lonE']
    theLatMin=myplot['latS']; theLatMax=myplot['latE']
    matrixToPlot = np.ma.masked_invalid(matrixToPlot)

# Set up the plot
    fig=pl.figure(num=None, figsize=(4,3), dpi=200, facecolor='w', edgecolor='k')
    ax1= fig.add_subplot(111)
#   map = Basemap(projection='ortho',lat_0=38,lon_0=14)
    map = Basemap(projection='merc',lat_0=38,lon_0=14,\
                                     llcrnrlon = -5.3, \
                                     llcrnrlat = 28.0, \
                                     urcrnrlon = 37, \
                                     urcrnrlat = 46.0, \
                                     resolution='l')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
#   map.drawcountries(linewidth=0.25)


# do the plot
    cs=map.pcolormesh(Lon,Lat,matrixToPlot,cmap=cmap,latlon='true')

# Set limits of scalar figure plotted
    theMin=myplot['vmin']
    theMax=myplot['vmax']
    pl.clim(theMin,theMax)

# Figure title
    theTitle = myplot['title']
    ax1.set_title(theTitle,fontsize=8)

#Insert annotations
    theUnit        = r'Units: '       + myplot['unit_label']              ; pl.annotate(theUnit ,xy=(0.1,0.35), xycoords='axes fraction' , fontsize=4);
    theStart       = r'From: '        + myplot['start']                   ; pl.annotate(theStart,xy=(0.1,0.30), xycoords='axes fraction',  fontsize=4);
    theEnd         = r'To__: '        + myplot['end']                     ; pl.annotate(theEnd ,xy=(0.1,0.25)  , xycoords='axes fraction', fontsize=4);
    theDepth       = r'Depth: '       + myplot['depth']                   ; pl.annotate(theDepth      ,xy=(0.1,0.20),xycoords='axes fraction',fontsize=4);

#   theAggregation = r'Aggregation: ' + myplot['aggregation']                 ; pl.annotate(theAggregation,xy=(0.1,0.20),xycoords='axes fraction', fontsize=4);
#   theDepth       = r'Depth: '       + str(myplot['depth']) + ' m '          ; pl.annotate(theDepth      ,xy=(0.1,0.15),xycoords='axes fraction',fontsize=4);

#-------------------------------------------------
    cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
    cbar.ax.tick_params(labelsize=5)
    return fig

def mapplot(map_dict, maskobj, fig, ax, ncolors=256, colormap='viridis'):
    """
    map_dict={'data': map2d, 'clim': [0,1], 'title': mystring }
    """
    background_color=(.9, .9, .9)
    cmap=pl.get_cmap(colormap,ncolors)

#    font=FontProperties()
    font_prop   = font_manager.FontProperties(fname='TitilliumWeb-Bold.ttf', size='xx-large', weight='bold')
    font_prop13 = font_manager.FontProperties(fname='TitilliumWeb-Regular.ttf', size=13)
    # font.set_name('TitilliumWeb')
    if (fig is None) or (ax is None):
        fig , ax = pl.subplots()
        fig.set_size_inches(10,10*map_obj.ymax/map_obj.xmax * (.85/.85))

    else:
        fig.clf()
        fig.add_axes(ax)

    ax.set_position([.10, .10, .85, .85])
    R1 = Rectangle(28.90338, 34.88078, 37.28154, 39.795254)#anatolia
    R2 = Rectangle(7.86973,11.0276,45.1554,46 )#nord italia
    R3 = Rectangle(11.7043,12.3246,42.3872,43.5016)# centro italia
    R4 = Rectangle(-0.532455,0.369799,40.9977,41.8014)# catalunya
    R5 = Rectangle(-6,-3.57754,39.011,42.0531)# centro spagna
    R6 = Rectangle(20.332,21.6854,40.3991,41.4642) #albania
    R7 = Rectangle(35.1627,35.9522,30.7444,33.0422)# israel
    R8 = Rectangle(28.2831,28.9034,45.1951,45.7487)# moldavia
    R9 = Rectangle(19.0258,19.588,42.0524,42.3803)# Skadarko Jezero
    R10 = Rectangle(27.7687,28.845,40.0348,40.2938)# 3 laghi vicino bursa, turchia
    R11 = Rectangle(29.2417,30.0876,39.9683,40.6135)# terzo lago
    R12 = Rectangle(31.7135,32.3184,31.0228,31.5918)#nilo
    RECTANGLE_LIST=[R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11] # R12]
    POLY=map_obj.coastpolygons
    for polygon in POLY:
        x,y=polygon
        ax.fill(x,y,color=background_color)
        if not is_in_boxes(x[0], y[0], RECTANGLE_LIST, map_obj):
            ax.plot(x,y,linewidth=0.7, color="0.4")


    vmin, vmax = map_dict['clim']
    map2d=map_dict['data']
    Zm = np.ma.masked_invalid(map2d)
    cs=map_obj.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm,cmap=cmap,latlon='true',vmin=vmin,vmax=vmax, shading="flat", ax=ax)
    parallels = np.arange(32.,46.,4)
    h_dict=map_obj.drawparallels(parallels,labels=[1,0,0,1],linewidth=0.5, fontsize=13, dashes=[1,2], ax=ax)
    #set_font(h_dict, font_prop13)
    h_dict = map_obj.drawparallels(parallels,labels=[1,0,0,1],fontsize=13, dashes=[6,900],ax=ax)
    set_invisible(h_dict)
    # draw meridians
    meridians = np.arange(-4.,40,4.)
    set_font(map_obj.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.5,fontsize=13, dashes=[1,2],ax=ax), font_prop13 ) #dashes=[6,900])
    set_invisible(map_obj.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, dashes=[6,900],ax=ax))

    nticks = 6
    cbar_ticks_list = np.linspace(vmin, vmax, nticks).tolist()
    cbar_ticks_labels = ["%g" % (t,)  for t in cbar_ticks_list ]


    cbar = map_obj.colorbar(cs,location='right',size="5%",pad="2%", ticks=cbar_ticks_list)
    cbar.ax.set_yticklabels(cbar_ticks_labels, fontproperties=font_prop13)


    t=ax.set_xlabel('Longitude',size= 'xx-large',fontweight='bold',labelpad = 20)#).set_fontsize(15)
    t.set_font_properties(font_prop)
    t=ax.set_ylabel('Latitude',size= 'xx-large',fontweight='bold',labelpad = 40)#).set_fontsize(15)
    t.set_font_properties(font_prop)

    ax.set_position([.10, .10, .85, .85])
    #title_1 = "%s, %s "  % (map_dict['longname'],map_dict['units'])
    #title_2 = "%s" %  map_dict['layer'].__repr__()
    #title_3 = "%s" % (map_dict['date']).strftime('%d - %m - %Y')
    #t1=ax.annotate(title_1 ,xy=(0.00,1.02), xycoords='axes fraction' , horizontalalignment='left')
    t2=ax.annotate(map_dict['title'] ,xy=(0.50,1.02), xycoords='axes fraction' , horizontalalignment='center')
    #t3=ax.annotate(title_3 ,xy=(1.00,1.02), xycoords='axes fraction' , horizontalalignment='right')
    #t1.set_font_properties(font_prop13)
    t2.set_font_properties(font_prop13)
    #t3.set_font_properties(font_prop13)


    return fig, ax