import pylab as pl

def mapplot(map_dict):
    clim = map_dict['clim']
    pl.imshow(map_dict['data'])
    pl.clim(clim[0], clim[1])
    pl.colorbar()
    pl.gca().invert_yaxis()
    pl.show()

