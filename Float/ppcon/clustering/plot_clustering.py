import numpy as np
#from dadapy.data import Data
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from dadapy.plot import plot_SLAn, plot_MDS, plot_matrix, get_dendrogram, plot_DecGraph


def plot_density_points(X, log_den):
    f, [ax1, ax2] = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'hspace': 0.05, 'wspace': 0})
    ax1.yaxis.set_major_locator(plt.NullLocator())
    ax1.xaxis.set_major_locator(plt.NullLocator())
    ax1.set_title('Estimated log densities')

    ax1.scatter(X[:, 0], X[:, 1], s=15., alpha=0.9, c=log_den, linewidths=0.0)
    ax2.yaxis.set_major_locator(plt.NullLocator())
    ax2.xaxis.set_major_locator(plt.NullLocator())
    ax2.set_title('Estimated log densities interpolated')
    ax2.tricontour(X[:, 0], X[:, 1], log_den, levels=10, linewidths=0.5, colors='k')
    fig2 = ax2.tricontourf(X[:, 0], X[:, 1], log_den, levels=250, alpha=0.9)

    plt.colorbar(fig2)
    plt.show()
    plt.close()
    return


def plot_adp_clustering(data):
    Nclus_m = len(data.cluster_centers)
    cmap = plt.get_cmap('gist_rainbow', Nclus_m)
    f, ax = plt.subplots(1, 1, figsize=(13, 10))
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.set_title('DPA assignation with halo')
    xdtmp = []
    ydtmp = []
    ldtmp = []
    xntmp = []
    yntmp = []
    for j in range(len(data.cluster_assignment)):
        if (data.cluster_assignment[j] != -1):
            xdtmp.append(data.X[j, 0])
            ydtmp.append(data.X[j, 1])
            ldtmp.append(data.cluster_assignment[j])
        else:
            xntmp.append(data.X[j, 0])
            yntmp.append(data.X[j, 1])

    plt.scatter(xdtmp, ydtmp, s=15., alpha=1.0, c=ldtmp, linewidths=0.0, cmap=cmap)
    plt.colorbar(ticks=range(Nclus_m))
    plt.clim(-0.5, Nclus_m - 0.5)
    plt.scatter(xntmp, yntmp, s=10., alpha=0.5, c='black', linewidths=0.0)
    plt.show()


def plot_float_coordinates(lat_values, lon_values, lat_limits_l=30, lat_limits_u=47, lon_limits_l=-3, lon_limits_u=37):
    lon_limits = [lon_limits_l, lon_limits_u]
    lat_limits = [lat_limits_l, lat_limits_u]
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


def plot_clustering_coordinates(lat_values, lon_values, clusters, lat_limits_l=30, lat_limits_u=47, lon_limits_l=-3,
                                lon_limits_u=37):
    lon_limits = [lon_limits_l, lon_limits_u]
    lat_limits = [lat_limits_l, lat_limits_u]
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
                              lat_values,
                              c=clusters)
    plt.ylim(lat_limits)
    plt.xlim(lon_limits)

    plt.show()
    plt.close()
