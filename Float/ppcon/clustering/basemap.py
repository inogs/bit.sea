import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd

lat_limits = [30, 47]
lon_limits = [-3, 37]


def plot_float_coordinates(lat_values, lon_values):
    fig = plt.figure(figsize=(8, 8))
    mediterranean_map = Basemap(llcrnrlon=lon_limits[0],
                                llcrnrlat=lat_limits[0],
                                urcrnrlon=lon_limits[1],
                                urcrnrlat=lat_limits[1],
                                resolution='h')
    mediterranean_map.drawmapboundary(fill_color='aqua')
    mediterranean_map.fillcontinents(color='coral', lake_color='aqua')
    mediterranean_map.drawcoastlines()

    samples_map = plt.scatter(lon_values,
                              lat_values)
    plt.ylim(lat_limits)
    plt.xlim(lon_limits)

    plt.show()
    plt.close()


plot_float_coordinates([33, 35], [-1, 23])
