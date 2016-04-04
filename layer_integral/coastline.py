import numpy as np

def get():
    coastline=np.load('Coastline.npy')
    c_lon=coastline['Lon']
    c_lat=coastline['Lat']
    return c_lon, c_lat