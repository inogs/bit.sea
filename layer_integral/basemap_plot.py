from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import numpy as np
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