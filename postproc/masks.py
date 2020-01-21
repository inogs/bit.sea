import numpy as np
class mesh():
    def __init__(self,jpi,jpj):
        self.jpi = jpi
        self.jpj = jpj
        self.lon = 0
        self.lat = 0
eps = 0.01; 
SatOrigMesh = mesh(733,253)
SatOrigMesh.lon = np.arange(-9.50, 36.25+eps, 1.0/16)
SatOrigMesh.lat = np.arange(30.25, 46.0 +eps, 1.0/16)

V4mesh = mesh(722,253)
V4mesh.lon = np.arange(-8.8125, 36.25+eps, 1.0/16)
V4mesh.lat = np.arange(30.1875, 45.9375+eps,1.0/16)

V1mesh = mesh(362,128)
V1mesh.lon = np.arange(-8.78125, 36.34375+eps, 1.0/8)
V1mesh.lat = np.arange(30.15620,46.0312+eps, 1.0/8)

KD490mesh = mesh(3308,1580)
KD490mesh.lon = np.linspace(-6.0, 36.500481, 3308, endpoint=True)
KD490mesh.lat = np.linspace(30.0, 45.998547, 1580, endpoint=True)

SAT1km_mesh = mesh(3308,1580)
SAT1km_mesh.lon = np.linspace(-6.0, 36.500481, 3308, endpoint=True)
SAT1km_mesh.lat = np.linspace(30.0, 45.998547, 1580, endpoint=True)

Mesh24 = mesh(1085,380)
Mesh24.lon = np.linspace(-8.875,36.291668,1085,endpoint=True)
Mesh24.lat = np.linspace(30.1875,45.979168,380,endpoint=True)

Mesh4 = mesh(182,65)
Mesh4.lon = np.linspace(-8.75  ,36.5  , 182,endpoint=True)
Mesh4.lat = np.linspace(30.1875,46.1875,  65,endpoint=True)

MeshOCCCI = mesh(1008, 384)
MeshOCCCI.lon = np.linspace(-5.97916650772095 , 35.9791679382324 , 1008,endpoint=True) #0.04166650772094993
MeshOCCCI.lat = np.linspace(30.0208339691162, 45.9791679382324,  384,endpoint=True) #0.04166603088379972
#approx values, because they have irregulare spacing
#max error about 1.e-4 deg in longitude, 3.e-4 in latitude