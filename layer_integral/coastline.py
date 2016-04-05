import numpy as np
import os

def get():
    localdir=os.path.dirname(os.path.abspath(__file__))
    coastline=np.load(localdir + '/Coastline.npy')
    c_lon=coastline['Lon']
    c_lat=coastline['Lat']
    return c_lon, c_lat