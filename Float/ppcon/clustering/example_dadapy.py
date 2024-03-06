import numpy as np
from dadapy.data import Data
import matplotlib.pyplot as plt
from dadapy.plot import plot_SLAn, plot_MDS, plot_matrix, get_dendrogram, plot_DecGraph


# Generate a simple 3D gaussian dataset
X1 = np.random.normal(0, 1, (1000, 12))
X2 = 3 + np.random.normal(0, 1, (1000, 12))

print(X1)

# print(type(X1))
X = np.concatenate((X1, X2))

# initialise the "Data" class with a set of coordinates
data = Data(X)

# compute distances up to the 100th nearest neighbour
data.compute_distances(maxk=100)

# compute the intrinsic dimension using the 2NN estimator
intrinsic_dim, _, intrinsic_dim_err = data.compute_id_2NN()

# check the value of the intrinsic dimension found
print(data.intrinsic_dim)

# compute the density of all points using a simple kNN estimator
log_den, log_den_err = data.compute_density_kNN(k=15)

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

# as an alternative, compute the density using a more sophisticated estimator
log_den, log_den_err = data.compute_density_PAk()


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

plt.hist(data.log_den)
plt.show()
plt.close()

# Compute the so-called decison graph and plot it (Note that the density has been computed in previous steps).
data.compute_DecGraph()
plot_DecGraph(data)


# find the statistically significant peaks of the density profile computed previously
data.compute_clustering_ADP(Z=1.5)

print(data.N_clusters)

Nclus_m=len(data.cluster_centers)
cmap = plt.get_cmap('gist_rainbow', Nclus_m)
f, ax = plt.subplots(1, 1, figsize = (13, 10))
ax.yaxis.set_major_locator(plt.NullLocator())
ax.xaxis.set_major_locator(plt.NullLocator())
ax.set_title('DPA assignation with halo')
xdtmp=[]
ydtmp=[]
ldtmp=[]
xntmp=[]
yntmp=[]
for j in range(len(data.cluster_assignment)):
    if (data.cluster_assignment[j]!=-1):
        xdtmp.append(data.X[j,0])
        ydtmp.append(data.X[j,1])
        ldtmp.append(data.cluster_assignment[j])
    else:
        xntmp.append(data.X[j,0])
        yntmp.append(data.X[j,1])

plt.scatter(xdtmp,ydtmp,s=15.,alpha=1.0, c=ldtmp,linewidths=0.0,cmap=cmap)
plt.colorbar(ticks=range(Nclus_m))
plt.clim(-0.5, Nclus_m-0.5)
plt.scatter(xntmp,yntmp,s=10.,alpha=0.5, c='black',linewidths=0.0)
plt.show()
