from basins.region import Polygon
from basins.basin import SimplePolygonalBasin, ComposedBasin
from basins.basin import Basin

DATADIR='/g100_scratch/userexternal/vdibiagi/articleBC_2023/SUPERCOASTAL/TXT/WKTcoords_12nm/'

filename = DATADIR + 'GolfoTrieste.txt'
LON, LAT = Basin.read_WKT_coords(filename)
got = Polygon(LON,LAT)
got = SimplePolygonalBasin('got', got, 'Gulf of Trieste 12nm')

filename = DATADIR + 'LagunaGradoMarano.txt'
LON, LAT = Basin.read_WKT_coords(filename)
lgm = Polygon(LON,LAT)
lgm = SimplePolygonalBasin('lgm', lgm, 'Grado-Marano lagoon 12nm')

filename = DATADIR + 'VenetoEst.txt'
LON, LAT = Basin.read_WKT_coords(filename)
vne = Polygon(LON,LAT)
vne = SimplePolygonalBasin('vne', vne, 'Eastern Veneto 12 nm')

filename = DATADIR + 'LagunaVenezia.txt'
LON, LAT = Basin.read_WKT_coords(filename)
lve = Polygon(LON,LAT)
lve = SimplePolygonalBasin('lve', lve, 'Venice lagoon 12nm')

filename = DATADIR + 'Po.txt'
LON, LAT = Basin.read_WKT_coords(filename)
poa = Polygon(LON,LAT)
poa = SimplePolygonalBasin('poa', poa, 'Po area 12nm')

filename = DATADIR + 'EmiliaRomagna.txt'
LON, LAT = Basin.read_WKT_coords(filename)
emr = Polygon(LON,LAT)
emr = SimplePolygonalBasin('emr', emr, 'Emilia-Romagna 12nm')

filename = DATADIR + 'Marche.txt'
LON, LAT = Basin.read_WKT_coords(filename)
mar = Polygon(LON,LAT)
mar = SimplePolygonalBasin('mar', mar, 'Marche 12 nm')

filename = DATADIR + 'Abruzzo.txt'
LON, LAT = Basin.read_WKT_coords(filename)
abr = Polygon(LON,LAT)
abr = SimplePolygonalBasin('abr', abr, 'Abruzzo 12nm')

filename = DATADIR + 'Molise.txt'
LON, LAT = Basin.read_WKT_coords(filename)
mol = Polygon(LON,LAT)
mol = SimplePolygonalBasin('mol', mol, 'Molise 12nm')

# add other polygons


# Composed basins

fvg = ComposedBasin('fvg',[got,lgm],'Friuli-Venezia Giulia 12 nm')
ven = ComposedBasin('ven',[vne,lve,poa],'Veneto 12 nm')

NAd = ComposedBasin('NAd',[got,lgm,vne,lve,poa,emr,mar],'NAdri 12 nm')


#SUBlist = NAd.basin_list
#for iSub, sub in enumerate(SUBlist):
#        print(iSub,sub.name)
