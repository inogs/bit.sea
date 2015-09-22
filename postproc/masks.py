import numpy as np
class mesh():
    def __init__(self,jpi,jpj):
        self.jpi = jpi
        self.jpj = jpj
        self.lon = 0
        self.lat = 0

SatOrigMesh = mesh(733,253)
SatOrigMesh.lon = np.arange(-9.50, 36.25, 1.0/16)
SatOrigMesh.lat = np.arange(30.25, 46.0 , 1.0/16)

V4mesh = mesh(722,253)
V4mesh.lon = np.arange(-8.8125, 36.25, 1.0/16)
V4mesh.lat = np.arange(30.1875, 45.9375,1.0/16)

V1mesh = mesh(362,128)
V1mesh.lon = np.arange(-8.78125, 36.34375, 1.0/8)
V1mesh.lat = np.arange(30.15620,46.0312, 1.0/8)