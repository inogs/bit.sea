from basins.region import Rectangle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import numpy as np
import matplotlib.dates as mdates
xlim=[-6.5,36.5]
ylim=[30,46]
xC=(xlim[0]+xlim[1])/2
yC=(ylim[0]+ylim[1])/2
font_s =  13
font_s2 = 3
label_s = 15
mapobj = Basemap(projection='merc',lat_0=xC,lon_0=yC,\
                                  llcrnrlon = xlim[0], \
                                  llcrnrlat = ylim[0], \
                                  urcrnrlon = xlim[1], \
                                  urcrnrlat = ylim[1], \
                                  area_thresh=None, \
                                  resolution='l')
background_color=(.9, .9, .9)
POLY=mapobj.coastpolygons
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

def is_in_boxes(x,y,Rectangle_list, map_obj):
    for rect in Rectangle_list:
        lonmin,latmin =map_obj(rect.lonmin, rect.latmin)
        lonmax,latmax =map_obj(rect.lonmax, rect.latmax)
        rect_meters = Rectangle(lonmin, lonmax, latmin, latmax)
        if rect_meters.is_inside(x,y):
            return True
    return False

def figure_generator(Profilelist):
    nP=len(Profilelist)
    Lon  = np.zeros((nP,), np.float32)
    Lat  = np.zeros((nP,), np.float32)    
    LonN = np.zeros((nP,), np.float32)
    LatN = np.zeros((nP,), np.float32)
    LonZ = np.zeros((nP,), np.float32)
    LatZ = np.zeros((nP,), np.float32)   
    for ip, p in enumerate(Profilelist):
        Lon[ip] = p.lon
        Lat[ip] = p.lat
        LonN[ip], LatN[ip] = mapobj(p.lon,p.lat)

    fig = pl.figure()
    fig.set_size_inches(12,15)
    fig.set_dpi(300)       
        # AXES FOR MAPS
    ax1 = pl.subplot2grid((7, 3), (0, 0), colspan=2)
    ax2 = pl.subplot2grid((7, 3), (0, 2))
    ax3 = pl.subplot2grid((7, 3), (1, 0), colspan=3)
    ax4 = pl.subplot2grid((7, 3), (2, 0), colspan=3)
    ax5 = pl.subplot2grid((7, 3), (3, 0), colspan=3)
    ax6 = pl.subplot2grid((7, 3), (4, 0), colspan=3)
    ax7 = pl.subplot2grid((7, 3), (5, 0), colspan=3)
#    ax8 = pl.subplot2grid((7, 3), (6, 0), colspan=3)
#    axes=[ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    axes=[ax1,ax2, ax3, ax4, ax5, ax6, ax7]
    
    for ax in axes: # abbassiamo un po
        l,b,w,t=ax.get_position().extents
        ax.set_position([l,b-0.05,w-l,t-b])

        
    
    for polygon in POLY:
        x,y=polygon
        ax1.fill(x,y,color=background_color)
        if not is_in_boxes(x[0], y[0], RECTANGLE_LIST, mapobj):
            ax1.plot(x,y,linewidth=0.7, color="0.4")

    l,b,w,t=ax1.get_position().extents
    ax1.set_position([l, b, 0.45, 0.45*mapobj.ymax/mapobj.xmax ])
    ax1.plot(LonN,LatN,'r.',markersize=font_s2)
    ax1.plot(LonN[0],LatN[0],'bo',markersize=font_s2)
    ax1.plot(LonN[0],LatN[0],'bx')    
    parallels = np.arange(32.,46.,4)
    mapobj.drawparallels(parallels,labels=[1,0,0,1],linewidth=0.5, fontsize=13, dashes=[1,2], ax=ax1)
    meridians = np.arange(-4.,40,8.)
    mapobj.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.5,fontsize=13, dashes=[1,2],ax=ax1)



    xmin=min(Lon)-2
    xmax=max(Lon)+2
    ymin=min(Lat)-2
    ymax=max(Lat)+2
    xC=(xmin+xmax)/2
    yC=(ymin+ymax)/2
    mapobj_zoom = Basemap(projection='merc',lat_0=yC,lon_0=xC,\
                                  llcrnrlon = xmin, \
                                  llcrnrlat = ymin, \
                                  urcrnrlon = xmax, \
                                  urcrnrlat = ymax, \
                                  area_thresh=None, \
                                  resolution='l')
    for ip, p in enumerate(Profilelist):
        LonZ[ip], LatZ[ip] = mapobj_zoom(p.lon,p.lat)

    POLY_zoom=mapobj_zoom.coastpolygons
    for polygon in POLY_zoom:
        x,y=polygon
        ax2.fill(x,y,color=background_color)
        if not is_in_boxes(x[0], y[0], RECTANGLE_LIST, mapobj):
            ax2.plot(x,y,linewidth=0.7, color="0.4")

    l,b,w,t=ax1.get_position().extents
    deltax=(1-b)/(mapobj_zoom.ymax/mapobj_zoom.xmax)
    if deltax>0.27 : deltax=0.27

    ax2.plot(LonZ,LatZ,'r.',markersize=font_s2)
    ax2.plot(LonZ[0],LatZ[0],'bo',markersize=font_s2)
    ax2.plot(LonZ[0],LatZ[0],'bx',markersize=font_s)

    ax2.set_position([0.90-deltax, b,deltax, deltax*mapobj_zoom.ymax/mapobj_zoom.xmax ])

    parallels = np.arange(30.,46.,2)
    h_dict=mapobj_zoom.drawparallels(parallels,labels=[1,0,0,1],linewidth=0.5, fontsize=13, dashes=[1,2], ax=ax2)
    meridians = np.arange(-4.,40,2.)
    mapobj_zoom.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.5,fontsize=13, dashes=[1,2],ax=ax2)
    ax2.set_position([0.90-deltax, b,deltax, deltax*mapobj_zoom.ymax/mapobj_zoom.xmax ])

    ax3.set_ylabel("FLOAT \n depth $[m]$",color = 'k', fontsize=font_s)
    ax4.set_ylabel("MODEL \ndepth $[m]$",color = 'k',fontsize=font_s)

    ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
    ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11])) 
    ax5.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
    ax6.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
    ax7.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,3,5,7,9,11]))
    ax7.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))    

    return fig, axes
