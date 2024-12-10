import numpy as np
import os, sys
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.basins import V2 as OGS


annaCoast = False

maskfile    = os.getenv("MASKFILE");
if maskfile is None :
    print("Error: Environment variable MASKFILE must be defined ")
    sys.exit(1)

if annaCoast:
    kcoastfile  = os.getenv( "KCOASTFILE");
    if kcoastfile is None:
        print("Error: Environment variable KCOASTFILE must be defined ")
        sys.exit(1)

TheMask = Mask.from_file(maskfile, ylevels_var_name="gphit", xlevels_var_name="glamt")
jpk,jpj,jpi = TheMask.shape


tmask   =  TheMask.mask
nav_lev =  TheMask.zlevels
Lon     =  TheMask.xlevels
Lat     =  TheMask.ylevels
area    =  TheMask.area
e3t     =  TheMask.dz


MEDmask      = tmask.copy()
MEDmask[:, Lon < -5.3] = False# Atlantic buffer 
if (nav_lev.max() > 200.0):
    mask200_2D   = TheMask.mask_at_level(200.0)
    #mask200_2D[Lon < -5.3] = False # Atlantic buffer
    mask200_3D = np.zeros((jpk,jpj,jpi),dtype=bool)
    for i in range(jpk):
        mask200_3D[i,:,:]=mask200_2D

if annaCoast:
    # Read k means class for coastal areas
    kcoastmask = np.load(kcoastfile)
    kmask1_2D =  np.zeros((jpj,jpi),dtype=bool)
    kmask1_2D[kcoastmask==1] = True
    kmask2_2D =  np.zeros((jpj,jpi),dtype=bool)
    kmask2_2D[kcoastmask==2] = True
    kmask1 =  np.zeros((jpk,jpj,jpi),dtype=bool)
    kmask2 =  np.zeros((jpk,jpj,jpi),dtype=bool)
    for i in range(jpk):
        kmask1[i,:,:] = kmask1_2D
        kmask2[i,:,:] = kmask2_2D

# Coastness must be defined one by one
COASTNESS_LIST=['coast', 'open_sea','everywhere']
if annaCoast: COASTNESS_LIST=['coast','open_sea','everywhere','coast1','coast2']
    
COASTNESS = np.ones((jpk,jpj,jpi) ,dtype=[(coast,bool) for coast in COASTNESS_LIST])
if (nav_lev.max() > 200.0):
    COASTNESS['open_sea']  =  mask200_3D;
    COASTNESS['coast']     = ~mask200_3D;

if annaCoast:
    COASTNESS['coast1']  =  kmask1;
    COASTNESS['coast2']  =  kmask2;

# Depth are defined by their names and bottom levels
DEPTHlist   =['shallow','mid','deep']
Bottom_list =[200, 600, 6000]

if annaCoast:
    DEPTHlist   =['shallow','deep']
    Bottom_list =[200, 6000]


DEPTH  = np.zeros((jpk,jpj,jpi),dtype=[(depth,bool) for depth in DEPTHlist])
tk_top = 0
for idepth, depth in enumerate(DEPTHlist):
    tk_bottom = TheMask.getDepthIndex(Bottom_list[idepth])+1
    DEPTH[depth][tk_top:tk_bottom ,:,:] = True
    tk_top = tk_bottom


Volume=np.zeros((jpk,jpj,jpi),dtype=np.double)
dZ     =np.ones ((jpk,jpj,jpi),dtype=np.double)
for i in range(jpk):
    Volume[i,:,:] = area*e3t[i]
    dZ[i,:,:]     = e3t[i]



SUBlist=[ sub.name for sub in OGS.P.basin_list ]
def SUB(sub):
    ''' sub is a string'''
    index= SUBlist.index(sub)
    basin = OGS.P.basin_list[index]
    s = SubMask(basin, TheMask)
    return s.mask



# SUB[med] is not applied by aveScan.
#If we want it we hav to apply to each subbasin


def read_Positions_for_Pointprofiles(filename):

    dtIN  = np.dtype([('Name','U100'), ('Lon',np.float32), ('Lat',np.float32)])
    dtOUT = np.dtype([('Name','U100'), ('Lon',np.float32), ('Lat',np.float32), ('i',int), ('j',int)])
    MeasPoints=np.loadtxt(filename, dtype=dtIN, skiprows=1,ndmin=1,delimiter="\t")
    
    nMeas=MeasPoints.size
    MeasPoints_OUT=np.zeros((nMeas),dtype=dtOUT)
    MeasPoints_OUT['Name']= MeasPoints['Name']
    MeasPoints_OUT['Lon' ]= MeasPoints['Lon' ]
    MeasPoints_OUT['Lat' ]= MeasPoints['Lat' ]
    for k in range(nMeas):
        i, j = TheMask.convert_lon_lat_wetpoint_indices(MeasPoints['Lon'][k],MeasPoints['Lat'][k])
#        DIST = (Lon - MeasPoints['Lon'][k])**2 + (Lat - MeasPoints['Lat'][k])**2 
#        ind=np.nonzero(DIST==DIST.min())# tuple  
#        j = ind[0]
#        i = ind[1]
        MeasPoints_OUT['i'][k] = i
        MeasPoints_OUT['j'][k] = j
    
    return MeasPoints_OUT
