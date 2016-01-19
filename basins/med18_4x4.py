
from basins.region import Rectangle
from basins.basin import SimplePolygonalBasin, ComposedBasin


def sequence(start_lon, start_lat, degrees):
    BL_point = [start_lon, start_lat]
    TR_point = [start_lon + degrees, start_lat + degrees]
    rect = Rectangle(BL_point[0], TR_point[0], BL_point[1], TR_point[1])
    return rect




LAT=[34., 35., 36., 40., 36.0, 41., 34.0, 39., 42., 32.5, 36.5, 40.5, 33, 31.0, 35., 39., 32.0, 32.]
LON=[-5.5, -1., 03., 03., 07.5, 06., 11.5, 10., 12., 16.0, 17., 16.0, 20., 24.0, 24., 23., 28.0, 32.]

np = 18; 
BASINLIST=[]
for i in range(18):
    rect = sequence(LON[i], LAT[i], 4.0)
    stringa = "sub_%02d" %i
    basin = SimplePolygonalBasin(stringa, rect)
    BASINLIST.append(basin)
    

P = ComposedBasin('P',BASINLIST,'med184x4')