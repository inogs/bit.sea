import pylab as pl
import numpy as np


def figure_generator(p):
    fig, axs = pl.subplots(2,3, facecolor='w', edgecolor='k')
    hsize=10
    vsize=12
    fig.set_size_inches(hsize,vsize)
    #figtitle = " date="+p.time.strftime('%Y/%m/%d')+" float="+p.name()
    #fig.set_title(figtitle)
    fig.subplots_adjust(hspace = 0.3, wspace=0.3)
    axs = axs.ravel()
    
    ax = axs[0]
    coastline=np.load('Coastline.npy')
    ax.plot(coastline['Lon'],coastline['Lat'], color='#000000',linewidth=0.5)
    ax.plot(p.lon,p.lat,'ro')
    ax.set_xticks(np.arange(-6,36))
    ax.set_yticks(np.arange(-30,46))
    ax.set_xlim([p.lon -2, p.lon+2])
    ax.set_ylim([p.lat -2, p.lat+2])
    bbox=ax.get_position()
    
    deltax, _ =bbox.size
    new_deltay = deltax* hsize/vsize
    bottom = bbox.ymax - new_deltay
    ax.set_position([bbox.xmin, bottom, deltax, new_deltay])     
    
    
    for ax in axs[1:]:
        ax.set_ylim(0,400)
        ax.yaxis.grid()
    
    for ax in [axs[2], axs[4], axs[5]]:
        ax.set_yticklabels([])
    return fig,axs

