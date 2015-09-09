import scipy.io.netcdf as NC
import numpy as np
import os, sys

def getDepthIndex(nav_lev, lev_mask):
    jk_m = 0
    for jk in range(nav_lev.__len__()):
        if nav_lev[jk] < lev_mask:
            jk_m = jk
    return jk_m


submaskfile = os.getenv("SUBMASKFILE");
maskfile    = os.getenv(   "MASKFILE"); 

if submaskfile is None or maskfile is None: 
    print "Error: Environment variables MASKFILE and SUBMASKFILE must be defined "
    sys.exit(1)

M=NC.netcdf_file(maskfile,"r")
jpi=M.dimensions['x']
jpj=M.dimensions['y']
jpk=M.dimensions['z']

tmask   = (M.variables['tmask'].data[0,:,:,:]).astype(np.bool)
nav_lev =  M.variables['nav_lev'].data
Lon     =  M.variables['glamt'].data[0,0,:,:]
Lat     =  M.variables['gphit'].data[0,0,:,:]

area    = (M.variables['e1t'].data[0,0,:,:])*(M.variables['e2t'].data[0,0,:,:])
e3t     =  M.variables['e3t'].data[0,:,0,0]

tk_m     = getDepthIndex(nav_lev,200.0)
MEDmask      = tmask.copy()
MEDmask[:, Lon < -5.3] = False# Atlantic buffer 

mask200_2D   = tmask[tk_m,:,:].copy()
#mask200_2D[Lon < -5.3] = False # Atlantic buffer
mask200_3D = np.zeros((jpk,jpj,jpi),dtype=np.bool)
for i in range(jpk):
    mask200_3D[i,:,:]=mask200_2D

COASTNESS_LIST=['coast','open_sea','everywhere']
struct=[]
for coast in COASTNESS_LIST:
    struct.append((coast,np.bool))
    
COASTNESS = np.ones((jpk,jpj,jpi),dtype=struct)
COASTNESS['coast']     = ~mask200_3D;
COASTNESS['open_sea']  =  mask200_3D;
#COASTNESS['everywhere'] = True;
DEPTHlist      =['shallow','deep']
struct=[]
for depth in DEPTHlist:
    struct.append((depth,np.bool))

DEPTH  = np.zeros((jpk,jpj,jpi),dtype=struct)
tk_1     = getDepthIndex(nav_lev,200.0)+1
#tk_2     = getDepthIndex(nav_lev,600.0)+1
#DEPTH['surf'   ][0        ,:,:] = True
DEPTH['shallow'][0:tk_1   ,:,:] = True
#DEPTH['mid'    ][tk_1:tk_2,:,:] = True
DEPTH['deep'   ][tk_1:    ,:,:] = True




Volume=np.zeros((jpk,jpj,jpi),dtype=np.double)
dZ     =np.ones ((jpk,jpj,jpi),dtype=np.double)
for i in range(jpk):
    Volume[i,:,:] = area*e3t[i]
    dZ[i,:,:]     = e3t[i]

M.close()


sbIN=NC.netcdf_file(submaskfile,"r")
Poly=__import__(sbIN.SubBasindef)
SUBlist=Poly.P #SUBlist=['alb','swm','nwm', 'tyr', 'adn','ads','aeg','ion','lev']
struct=[]
for sub in SUBlist:
    struct.append((sub,np.bool))
struct.append(('med',np.bool))

SUB=np.array(np.zeros((jpk,jpj,jpi)),dtype=struct)

for sub in SUBlist:
    SUB[sub] = sbIN.variables[sub].data
sbIN.close()
SUB['med']=MEDmask

SUBlist.append('med')


def read_Positions_for_Pointprofiles(filename):

    dtIN  = np.dtype([('Name','S20'), ('Lon',np.float32), ('Lat',np.float32)])
    dtOUT = np.dtype([('Name','S20'), ('Lon',np.float32), ('Lat',np.float32), ('i',np.int), ('j',np.int)])  
    MeasPoints=np.loadtxt(filename, dtype=dtIN, skiprows=1,ndmin=1)
    
    nMeas=MeasPoints.size
    MeasPoints_OUT=np.zeros((nMeas),dtype=dtOUT)
    MeasPoints_OUT['Name']= MeasPoints['Name']
    MeasPoints_OUT['Lon' ]= MeasPoints['Lon' ]
    MeasPoints_OUT['Lat' ]= MeasPoints['Lat' ]
    for k in xrange(nMeas):
        DIST = (Lon - MeasPoints['Lon'][k])**2 + (Lat - MeasPoints['Lat'][k])**2 
        ind=np.nonzero(DIST==DIST.min())# tuple  
        j = ind[0]
        i = ind[1]
        MeasPoints_OUT['i'][k] = i
        MeasPoints_OUT['j'][k] = j
    
    return MeasPoints_OUT
